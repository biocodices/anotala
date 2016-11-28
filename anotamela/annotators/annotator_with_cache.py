import sys
import time
import logging
from concurrent.futures import ThreadPoolExecutor

from tqdm import tqdm

from anotamela.helpers import grouped
from anotamela.cache import RedisCache


logger = logging.getLogger(__name__)


class AnnotatorWithCache():
    """
    Abstract Base Class for Annotators like DbSNP and MyVariant.
    The point of this class is to provide a shared logic of caching the
    responses of remote APIs and of parallelizing requests.

    To use this class, create a new annotator class that implements
    `_query()` and `_key()` methods. Optionally, you can override
    `_batch_query()` for services that you want to parallelize in a different
    way or not parallelize at all.
    """
    def __init__(self, cache='redis'):
        self.name = self.__class__.__name__

        if cache == 'redis':
            self.cache = RedisCache()

    def annotate(self, ids, parallel=10, sleep_time=10, use_cache=True,
                 use_web=True):
        """
        Annotate one or more IDs. If use_cache is set, cached responses will
        be prioritized and web will be only used for the missing annotations.
        If use_web is set, web annotation will be used as a source, if not,
        you will only get the annotations for the cached IDs. If use_cache is
        False, then all the IDs will be fetched from web (this can be used to
        update the cached data).

        Some annotators implement extra parsing of the web response with
        parse_info_dict(). This is enabled by default, but you can disable it
        with parse_data=False.
        """
        if isinstance(ids, str) or isinstance(ids, int):
            ids = [ids]
        ids = set(ids)

        info_dict = {}

        if use_cache:
            # FIXME: Convert ids to KEYS before calling cache
            # keys = [self._key(key) for key in keys]
            info_dict.update(self._cache_get(ids))
            ids = ids - info_dict.keys()
        else:
            logger.info(' Not using cache')

        if use_web:
            if ids:
                info_from_api = self._batch_query(ids, parallel, sleep_time)
                info_dict.update(info_from_api)
                ids = ids - info_dict.keys()
        else:
            logger.info(' Not using web')

        if ids:
            logger.info(" No info for %s IDs" % len(ids))

        if parse_data:
            self.parse_info_dict(info_dict)

        return info_dict

    def _query_and_set_cache(self, id_):
        """Make a query for a single id, cache the response, and return it."""
        response = self._query(id_)
        if response:
            # FIXME: Convert ids to KEYS before calling cache
            self._cache_set({id_: response})
        return response

    def _query(self):
        raise NotImplementedError()

    def _key(self):
        raise NotImplementedError()

    def _batch_query(self, ids, parallel, sleep_time):
        """
        Query a group of IDs using <parallel> threads and cache the responses.
        It returns a dict with the queried info per ID.
        """
        grouped_ids = list(grouped(ids, parallel))

        msg = ('🌎 Get {} {} entries in {} batches ({} items/batch), '
               'sleeping {}s between batches')
        msg = msg.format(len(ids), self.name, len(grouped_ids), parallel,
                         sleep_time)
        logger.info(msg)

        with ThreadPoolExecutor(max_workers=parallel) as executor:
            sys.stdout.flush()  # Hack for tqdm progress bar correct display
            iterator = enumerate(tqdm(grouped_ids, total=len(grouped_ids)))
            for i, ids_group in iterator:
                if i > 0:
                    time.sleep(sleep_time)

                executor.map(self._query_and_set_cache, ids_group)

        # After caching all the responses, use the logic in #_cache_get()
        # to bring them from cache. The jsons will be json-loaded, the xml
        # will be left as they are, etc.

        # FIXME: Convert ids to KEYS before calling cache
        return self._cache_get(ids, verbose=False)

    def parse_info_dict(self, info_dict):
        # Override this for extra parsing logic of the cached raw data
        pass
