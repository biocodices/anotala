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


class NoProxiesException(Exception):
    pass


class ParallelAnnotator(AnnotatorWithCache):
    """
    Base class for annotators that have a _url(id_) method or alternatively
    a _query(id) method. This class implements a 'manual' parallelization with
    threads.

    - Modify the class variables BATCH_SIZE and SLEEP_TIME to tweak the
      parallelization behavior.
    - Add a _parse_annotation(raw) method to indicate how raw responses from
      the URLs should be parsed.

    Example:

        Class MyAnnotator(ParallelAnnotator):
            BACTH_SIZE = 500

            def _url(self, id_):
                return 'http://annotation-service.com/?q={}'.format(id_)

            @staticmethod
            def _parse_annotation(raw):
                return raw.get('some-key')

    """
    BATCH_SIZE = 10
    SLEEP_TIME = 2.5
    RANDOMIZE_SLEEP_TIME = False
    USER_AGENT_GENERATOR = UserAgent()

    @property
    def random_user_agent(self):
        return self.USER_AGENT_GENERATOR.random

    def _query(self, id_):
        url = self._url(id_)
        headers = {'User-Agent': self.random_user_agent}
        response = requests.get(url, headers=headers, proxies=self.proxies)
        if not response.ok:
            logger.warning('{} response: {}'.format(response.status_code, url))
            response.raise_for_status()

        if self.ANNOTATIONS_ARE_JSON:
            return response.json()
        else:
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

        if self.proxies:
            logger.info('{} using proxies: {}'.format(self.name, self.proxies))
        elif self.proxies == {}:
            logger.warning('{} not using proxies!'.format(self.name))
        elif self.proxies is None:  # i.e. different from emtpy dict {}
            message = ('No proxies set for {}. Please set self.proxies to '
                       'avoid a possible ban or explicitely set proxies as an '
                       'empty dict (self.proxies={{}}) if you want to proceed '
                       'anyway.'.format(self.name))
            raise NoProxiesException(message)

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

