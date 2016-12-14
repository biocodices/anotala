import sys
import time
import logging
from concurrent.futures import ThreadPoolExecutor

import requests
from numpy.random import random_sample
from fake_useragent import UserAgent
from tqdm import tqdm

from anotamela.annotators.base_classes import AnnotatorWithCache
from anotamela.helpers import grouped


logger = logging.getLogger(__name__)


class ParallelAnnotator(AnnotatorWithCache):
    """
    Base class for annotators that have a _query(id) method but no
    parallelization. This class implements a 'manual' parallelization with
    Threads. Modify class variables BATCH_SIZE and SLEEP_TIME to tweak the
    parallelization behavior. Set PROXIES dict for requests.get() if you need
    that.
    """
    BATCH_SIZE = 10
    SLEEP_TIME = 10
    PROXIES = {}
    RANDOMIZE_SLEEP_TIME = False

    def _query(self):
        raise NotImplementedError()

    def _query_with_random_user_agent(self, url):
        if not hasattr(self, 'user_agent_generator'):
            self.user_agent_generator = UserAgent()
        headers = {'user-agent': self.user_agent_generator.random}
        response = requests.get(url, headers=headers,
                                proxies=self.PROXIES or {})
        if not response.ok:
            logger.warning('{} response: {}'.format(response.status_code, url))
            response.raise_for_status()
        return response.text

    def _batch_query(self, ids):
        """
        Query a group of IDs using <class>.BATCH_SIZE threads, sleeping
        <class>.SLEEP_TIME between batches.

        Set <class>.RANDOMIZE_SLEEP_TIME = True to make the sleep time between
        batches random, in order to fool crawler detection.

        Yields tuples of (id_, annotation).
        """
        grouped_ids = list(grouped(ids, self.BATCH_SIZE))
        msg = ('{}: get {} entries in {} batches '
               '({} items/batch & sleep {}s between batches)')
        logger.info(msg.format(self.name, len(ids), len(grouped_ids),
                               self.BATCH_SIZE, self.SLEEP_TIME))

        if self.PROXIES:
            logger.info('{} using proxies: {}'.format(self.name, self.PROXIES))

        with ThreadPoolExecutor(max_workers=self.BATCH_SIZE) as executor:
            sys.stdout.flush()  # Hack to display tqdm progress bar correctly
            iterator = tqdm(grouped_ids, total=len(grouped_ids))
            for i, batch_of_ids in enumerate(iterator):
                if i > 0:
                    time.sleep(self.sleep_time)
                annotations = executor.map(self._query, batch_of_ids)
                yield dict(zip(batch_of_ids, annotations))

    @property
    def sleep_time(self):
        if self.RANDOMIZE_SLEEP_TIME:
            sleep_time = self.randomize(self.SLEEP_TIME)
        else:
            sleep_time = self.SLEEP_TIME

        return sleep_time

    @staticmethod
    def randomize(n):
        return n + random_sample() * n

