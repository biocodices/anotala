import sys
import time
import logging
from concurrent.futures import ThreadPoolExecutor

from tqdm import tqdm

from anotamela.annotators import AnnotatorWithCache
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

    def _query(self):
        raise NotImplementedError()

    def _batch_query_and_cache(self, ids):
        """
        Query a group of IDs using <class>.BATCH_SIZE threads, sleep
        <class>.SLEEP_TIME between batches, and cache the responses.

        Returns a dict with the queried info per ID.
        """
        grouped_ids = list(grouped(ids, self.BATCH_SIZE))
        msg = '{}: get {} entries in {} batches ({} items/batch)'
        logger.info(msg.format(self.name, len(ids), len(grouped_ids),
                               self.BATCH_SIZE))

        annotations = {}
        with ThreadPoolExecutor(max_workers=self.BATCH_SIZE) as executor:
            sys.stdout.flush()  # Hack to display tqdm progress bar correctly
            iterator = tqdm(grouped_ids, total=len(grouped_ids))
            for i, ids_group in enumerate(iterator):
                if i > 0:
                    time.sleep(self.SLEEP_TIME)
                group_annotations = executor.map(self._query, ids_group)

        group_annotations = {id_: annotation for id_, annotation
                             in zip(ids_group, group_annotations)
                             if annotation}
        self.cache.set(group_annotations, namespace=self.SOURCE_NAME)
        annotations.update(group_annotations)

        return annotations

