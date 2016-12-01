import sys
import time
import logging
from concurrent.futures import ThreadPoolExecutor

from tqdm import tqdm

from anotamela.helpers import grouped
from anotamela.cache import RedisCache, PostgresCache


logger = logging.getLogger(__name__)


class AnnotatorWithCache():
    """
    Abstract Base Class for Annotators with Cache. The point of this class is
    to provide a shared logic of caching the responses of remote APIs and of
    parallelizing requests.

    To use this class, create a new annotator class that implements a
    `_query()` method and has a SOURCE_NAME class variable. Optionally, you can
    override `_batch_query()` for services that you want to parallelize in a
    special way.

    SOURCE_NAME will work as a namespace or tablename according to each Cache
    used. For instance, cache='redis' will use the SOURCE_NAME as a prefix
    to the ID being cached (e.g. 'dbsnp:rs123'), whereas cache='postgres' will
    use the SOURCE_NAME as the name of the table where to store the info.
    """
    AVAILABLE_CACHES = {
            'redis': RedisCache,
            'postgres': PostgresCache
        }

    def __init__(self, cache='redis', **cache_kwargs):
        self.name = self.__class__.__name__

        try:
            cache_class = self.AVAILABLE_CACHES[cache]
        except KeyError as error:
            known_caches = ', '.join(self.AVAILABLE_CACHES.keys())
            msg = 'Unknown cache "{}". I only know: {}'.format(cache,
                                                               known_caches)
            raise ValueError(msg).with_traceback(error.__traceback__)

        self.cache = cache_class(**cache_kwargs)

    def annotate(self, ids, parallel=10, sleep_time=10, use_cache=True,
                 use_web=True, parse_data=True):
        """
        Annotate one or more IDs and return an info dictionary with one key per
        passed ID. If use_cache is set, cached responses will be prioritized
        and web will only be used for the missing annotations. If use_web is
        set, web annotation will be used as a source, if not, you will only get
        the annotations for the cached IDs. If use_cache is False, then all the
        IDs will be fetched from web (this can be used to update the cached
        data).

        Annotators can implement extra parsing of the data with
        parse_annotations_dict(). This parsing is enabled by default, but you
        can disable it with parse_data=False and get the raw web/cached responses.
        """
        ids = self._set_of_string_ids(ids)

        annotations_dict = {}

        if use_cache:
            cached_data = self.cache.get(ids, namespace=self.SOURCE_NAME)
            annotations_dict.update(cached_data)
            ids = ids - annotations_dict.keys()
        else:
            logger.info('Not using cache')

        if use_web:
            if ids:
                info_from_api = self._batch_query(ids, parallel, sleep_time)
                annotations_dict.update(info_from_api)
                ids = ids - annotations_dict.keys()
        else:
            logger.info('Not using web')

        if ids:
            logger.info('No info for {} IDs'.format(len(ids)))

        if parse_data:
            self.parse_annotations_dict(annotations_dict)

        return annotations_dict

    def _query_and_set_cache(self, id_):
        """Make a query for a single id, cache the response, and return it."""
        response = self._query(id_)
        if response:
            self.cache.set({id_: response}, namespace=self.SOURCE_NAME)
        return response

    def _query(self):
        raise NotImplementedError()

    def _batch_query(self, ids, parallel, sleep_time):
        """
        Query a group of IDs using <parallel> threads and cache the responses.
        It returns a dict with the queried info per ID.
        """
        grouped_ids = list(grouped(ids, parallel))

        msg = '{}: get {} entries in {} batches ({} items/batch)'
        logger.info(msg.format(self.name, len(ids), len(grouped_ids), parallel))

        annotations_dict = {}
        with ThreadPoolExecutor(max_workers=parallel) as executor:
            sys.stdout.flush()  # Hack to display tqdm progress bar correctly
            iterator = enumerate(tqdm(grouped_ids, total=len(grouped_ids)))
            for i, ids_group in iterator:
                if i > 0:
                    time.sleep(sleep_time)

                annotations = executor.map(self._query_and_set_cache, ids_group)
                annotations_dict.update(zip(ids_group, annotations))

        return annotations_dict

    def parse_annotations_dict(self, info_dict):
        # Override this for extra parsing logic
        pass

    @staticmethod
    def _set_of_string_ids(ids):
        if isinstance(ids, str) or isinstance(ids, int):
            ids = [ids]
        return set(str(id_) for id_ in set(ids))

