import logging
from concurrent.futures import (
        ProcessPoolExecutor,
        as_completed
    )

from anotamela.cache import Cache, create_cache


logger = logging.getLogger(__name__)


class AnnotatorWithCache():
    """
    Abstract Base Class for annotators. The point of this class is to provide
    a shared logic of annotating a list of IDs with the cache and/or with a web
    service, and then cache the new data retrieved from the latter.

    To use this class, create a new annotator class that has:

        - a SOURCE_NAME class variable that will work as a namespace or
          tablename according to each Cache used. For instance, RedisCache
          will use the SOURCE_NAME as a prefix to the ID being cached (e.g.
          'dbsnp:rs123'), whereas PostgresCache will use the SOURCE_NAME as the
          name of the table where to store the info.

        - either a `_query()` method to fetch a single ID's data --the
          annotation will be parallelized with multithread calls to that
          method-- or a `_batch_query_and_cache()` method to fetch a group of
          IDs, in case you want to implement parallelization in a different way

        - an optional @classmethod or @staticmethod _parse_annotation() that
          takes the annotation for one id and transforms it in any way.

        - an optional ANNOTATIONS_ARE_JSON class variable if the annotations
          retrieved from the web will be in JSON format. This lets
          PostgresCache know it can create a JSONB field in the database.

    """
    def __init__(self, cache='redis', **cache_kwargs):
        """
        Initialize with a cache name ('redis', 'postgres') or a Cache instance
        (RedisCache, PostgresCache). Extra kwargs can be passed to the cache
        initializer.
        """
        self.name = self.__class__.__name__

        if isinstance(cache, Cache):
            self.cache = cache
        else:
            self.cache = create_cache(cache, **cache_kwargs)

    def annotate(self, ids, use_cache=True, use_web=True, parse_data=True):
        """
        Annotate one or more IDs and return an info dictionary with the passed
        IDs as keys. If use_cache=True, cached responses will be prioritized
        and web will only be used for the missing annotations. If use_web=True,
        web annotation will be used as a source, if not, you will only get
        the annotations for the cached IDs. If use_cache=False, then all the
        IDs will be fetched from web (this can be used to update the cached
        data).

        Annotators can implement extra parsing of the data with
        _parse_annotation(). This parsing is enabled by default, but you can
        disable it with parse_data=False and get the raw responses.
        """
        ids = self._set_of_string_ids(ids)
        total_count = len(ids)
        msg = '{} annotating {} ids'.format(self.name, total_count)
        logger.info(msg)

        annotations = {}
        if use_cache:
            logger.info('Get info from cache')
            cached_data = self.cache.get(
                    ids,
                    namespace=self.SOURCE_NAME,
                    load_as_json=hasattr(self, 'ANNOTATIONS_ARE_JSON')
                )
            annotations.update(cached_data)
            ids = ids - annotations.keys()
        else:
            logger.info('Not using cache')

        if use_web:
            if ids:
                logger.info('Get info from web for {} IDs'.format(len(ids)))
                for batch_annotations in self._batch_query(ids):
                    self.cache.set(
                            batch_annotations,
                            namespace=self.SOURCE_NAME,
                            save_as_json=hasattr(self, 'ANNOTATIONS_ARE_JSON')
                        )
                    annotations.update(batch_annotations)
                    ids = ids - batch_annotations.keys()
        else:
            logger.info('Not using web')

        if ids:
            msg = '{} found info for {}/{} ({:.2%}) IDs'
            found_count = total_count - len(ids)
            logger.info(msg.format(self.__class__.__name__, found_count,
                                   total_count, found_count/total_count))

        if parse_data and hasattr(self, '_parse_annotation'):
            logger.info('Parsing {} annotations'.format(len(annotations)))
            annotations = self._parse_annotations(annotations)

        return annotations

    def _parse_annotations(self, annotations):
        """Parse a dict of annotations in parallel. Return a dictionary with
        the same keys and the parsed annotations."""
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

