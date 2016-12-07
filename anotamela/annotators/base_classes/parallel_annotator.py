import sys
import time
import logging
from concurrent.futures import ThreadPoolExecutor

from numpy.random import random_sample
from tqdm import tqdm

from anotamela.annotators.base_classes import AnnotatorWithCache
from anotamela.helpers import grouped


logger = logging.getLogger(__name__)


class ParallelAnnotator(AnnotatorWithCache):
    """
    Base class for annotators that have a _query(id) method but no
    parallelization. This class implements a 'manual' parallelization with
    Threads. Modify <class>.BATCH_SIZE and <class>.SLEEP_TIME to tweak the
    parallelization behavior.
    """
    BATCH_SIZE = 10
    SLEEP_TIME = 10
    RANDOMIZE_SLEEP_TIME = False

    def _query(self):
        raise NotImplementedError()

    def _batch_query_and_cache(self, ids):
        """
        Query a group of IDs using <class>.BATCH_SIZE threads, sleep
        <class>.SLEEP_TIME between batches, and cache the responses.

        Set <class>.RANDOMIZE_SLEEP_TIME = True to make the sleep time between
        batches random, in order to fool crawler detection.

        Returns a dict with the queried info per ID.
        """
        grouped_ids = list(grouped(ids, self.BATCH_SIZE))
        msg = ('{}: get {} entries in {} batches '
               '({} items/batch & sleep {} between batches)')
        logger.info(msg.format(self.name, len(ids), len(grouped_ids),
                               self.BATCH_SIZE, self.SLEEP_TIME))

        annotations = {}
        with ThreadPoolExecutor(max_workers=self.BATCH_SIZE) as executor:
            sys.stdout.flush()  # Hack to display tqdm progress bar correctly
            iterator = tqdm(grouped_ids, total=len(grouped_ids))
            for i, ids_group in enumerate(iterator):
                if i > 0:
                    time.sleep(self.sleep_time)
                group_annotations = executor.map(self._query, ids_group)
                group_annotations = {id_: annotation for id_, annotation
                                     in zip(ids_group, group_annotations)
                                     if annotation}
                self.cache.set(group_annotations, namespace=self.SOURCE_NAME)
                annotations.update(group_annotations)

        return annotations

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

