import sys
import time
import logging
from concurrent.futures import (
        ThreadPoolExecutor,
        ProcessPoolExecutor,
        as_completed
    )

from tqdm import tqdm

from anotamela.helpers import grouped
from anotamela.cache import Cache, RedisCache, PostgresCache


logger = logging.getLogger(__name__)


class AnnotatorWithCache():
    """
    Abstract Base Class for Annotators with Cache. The point of this class is
    to provide a shared logic of caching the responses of remote APIs and of
    parallelizing requests.

    To use this class, create a new annotator class that has:
        - SOURCE_NAME [and ANNOTATIONS_ARE_JSON] class variables
        - either a `_query()` method to fetch a single ID's data --the
          annotation will be parallelized with multithread calls to that
          method-- or a `_batch_query_and_cache()` method to fetch a group of
          IDs, in case you want to implement parallelization in a different way
        - an optional @classmethod or @staticmethod _parse_annotation() that
          takes the annotation for one id and transforms it in any way.

    SOURCE_NAME will work as a namespace or tablename according to each Cache
    used. For instance, cache='redis' will use the SOURCE_NAME as a prefix
    to the ID being cached (e.g. 'dbsnp:rs123'), whereas cache='postgres' will
    use the SOURCE_NAME as the name of the table where to store the info.

    ANNOTATIONS_ARE_JSON should be set only if the annotations will be JSON
    formatted. This lets PostgresCache know if the columns should be of type
    JSONB. RedisCache ignores this value.
    """
    AVAILABLE_CACHES = {
            'redis': RedisCache,
            'postgres': PostgresCache
        }

    def __init__(self, cache='redis', **cache_kwargs):
        """
        Initialize this class with a cache name ('redis', 'postgres') or a
        Cache instance (RedisCache, PostgresCache). Extra kwargs can be
        passed to the cache initializer.
        """
        self.name = self.__class__.__name__

        if isinstance(cache, Cache):
            self.cache = cache
        else:
            try:
                cache_class = self.AVAILABLE_CACHES[cache]
            except KeyError as error:
                known_caches = ', '.join(self.AVAILABLE_CACHES.keys())
                msg = 'Unknown cache "{}". I only know: {}'.format(cache,
                                                                   known_caches)
                raise ValueError(msg).with_traceback(error.__traceback__)
            else:
                self.cache = cache_class(**cache_kwargs)

        # Tell the cache if the annotations are going to be JSON
        # PostgresCache will know it should create JSONB fields in the db:
        self.cache.SAVE_AS_JSON = hasattr(self, 'ANNOTATIONS_ARE_JSON')

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
        _parse_annotation(). This parsing is enabled by default, but you can
        disable it with parse_data=False and get the raw web/cached responses.
        """
        ids = self._set_of_string_ids(ids)
        logger.info('{} annotating {} ids'.format(self.__class__.__name__,
                                                  len(ids)))

        annotations = {}

        if use_cache:
            cached_data = self.cache.get(ids, namespace=self.SOURCE_NAME)
            annotations.update(cached_data)
            ids = ids - annotations.keys()
        else:
            logger.info('Not using cache')

        if use_web:
            if ids:
                info_from_api = self._batch_query_and_cache(ids, parallel,
                                                            sleep_time)
                annotations.update(info_from_api)
                ids = ids - annotations.keys()
        else:
            logger.info('Not using web')

        if ids:
            msg = '{} had no info for {} IDs'
            logger.info(msg.format(self.__class__.__name__, len(ids)))

        if parse_data:
            annotations = self._parse_annotations(annotations)

        return annotations

    def _query(self):
        raise NotImplementedError()

    def _batch_query_and_cache(self, ids, parallel, sleep_time):
        """
        Query a group of IDs using <parallel> threads and cache the responses.
        It returns a dict with the queried info per ID.
        """
        grouped_ids = list(grouped(ids, parallel))

        msg = '{}: get {} entries in {} batches ({} items/batch)'
        logger.info(msg.format(self.name, len(ids), len(grouped_ids), parallel))

        annotations = {}
        with ThreadPoolExecutor(max_workers=parallel) as executor:
            sys.stdout.flush()  # Hack to display tqdm progress bar correctly
            iterator = enumerate(tqdm(grouped_ids, total=len(grouped_ids)))
            for i, ids_group in iterator:
                if i > 0:
                    time.sleep(sleep_time)
                group_annotations = executor.map(self._query, ids_group)

        group_annotations = {id_: annotation for id_, annotation
                                in zip(ids_group, group_annotations)
                                if annotation}
        self.cache.set(group_annotations, namespace=self.SOURCE_NAME)
        annotations.update(group_annotations)

        return annotations

    def _parse_annotations(self, annotations):
        """
        Parse a dict of annotations in parallel. Return a dictionary with the
        same keys and the parsed annotations.
        """
        if not hasattr(self, '_parse_annotation'):
            return annotations

        for id_, annotation in annotations.items():
            parsed_annotations = {}
            with ProcessPoolExecutor() as executor:
                future_to_id = {}
                for id_, annotation in annotations.items():
                    future = executor.submit(self._parse_annotation, annotation)
                    future_to_id[future] = id_
                for future in as_completed(future_to_id):
                    id_ = future_to_id[future]
                    parsed_annotations[id_] = future.result()

            return parsed_annotations

    @staticmethod
    def _set_of_string_ids(ids):
        if isinstance(ids, str) or isinstance(ids, int):
            ids = [ids]
        return set(str(id_) for id_ in set(ids))

