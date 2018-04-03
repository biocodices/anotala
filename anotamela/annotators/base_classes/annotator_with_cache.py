import logging

from anotamela.annotators.base_classes import Annotator
from anotamela.cache import Cache, create_cache


logger = logging.getLogger(__name__)


class AnnotatorWithCache(Annotator):
    """
    Abstract Base Class for web annotators with cache. The point of this class
    is to provide a shared logic of annotating a list of IDs with the cache
    and with a web service, and then cache the new data retrieved from the
    latter.

    To use this class, create a new annotator class that has:

        - a SOURCE_NAME class variable that will work as a namespace or
          tablename according to each Cache used. For instance, RedisCache
          will use the SOURCE_NAME as a prefix to the ID being cached (e.g.
          'dbsnp:rs123'), whereas PostgresCache will use the SOURCE_NAME as the
          name of the table where to store the info.

        - either a `_query()` method to fetch a single ID's data --the
          annotation will be parallelized with multithread calls to that
          method-- or a `_batch_query()` method to fetch a group of
          IDs, in case you want to implement parallelization in a different way

        - an optional @classmethod or @staticmethod _parse_annotation() that
          takes the annotation for one id and transforms it in any way.

        - an optional ANNOTATIONS_ARE_JSON=True class variable if the
          annotations retrieved from the web will be in JSON format. This lets
          PostgresCache know it can create a JSONB field in the database.

    """
    ANNOTATIONS_ARE_JSON = False  # Default to be overriden
    SOURCE_NAME = ''

    def __init__(self, cache='redis', proxies=None, **cache_kwargs):
        """
        Initialize with a cache name ('redis', 'postgres') or a Cache instance
        (RedisCache, PostgresCache). Extra kwargs can be passed to the cache
        initializer.

        Set proxies as a dict like {'http': 'socks5://localhost:9050'} or as
        an empty dict in case you're using a ParallelAnnotator.
        """
        self.name = self.__class__.__name__
        self.proxies = proxies
        self.cache_kwargs = cache_kwargs

        if isinstance(cache, Cache):
            self.cache = cache
        else:
            self.cache = create_cache(cache, **cache_kwargs)

    def __repr__(self):
        msg = "{}(cache={}, **{})"
        return msg.format(self.name, self.cache, self.cache_kwargs)

    def _annotate_many_ids(self, ids_to_annotate, use_cache=True, use_web=True):
        """
        Annotate one or more IDs and return an info dictionary with the passed
        IDs as keys. If use_cache=True, cached responses will be prioritized
        and web will only be used for the missing annotations. If use_web=True,
        web annotation will be used as a source, if not, you will only get
        the annotations for the cached IDs. If use_cache=False, then all the
        IDs will be fetched from web (this can be used to update the cached
        data).
        """
        annotations = {}
        if use_cache:
            logger.info('{} get info from cache'.format(self.name))
            cached_data = self.cache.get(
                    ids_to_annotate,
                    namespace=self.SOURCE_NAME,
                    as_json=self.ANNOTATIONS_ARE_JSON
                )
            annotations.update(cached_data)
            ids_to_annotate = ids_to_annotate - annotations.keys()
        else:
            logger.info('{} not using cache'.format(self.name))

        if use_web:
            if ids_to_annotate:
                logger.info('{} get info from web for {} IDs'
                            .format(self.name, len(ids_to_annotate)))
                for batch_annotations in self._batch_query(ids_to_annotate):
                    self.cache.set(
                            batch_annotations,
                            namespace=self.SOURCE_NAME,
                            as_json=self.ANNOTATIONS_ARE_JSON
                        )
                    annotations.update(batch_annotations)
                    ids_to_annotate = ids_to_annotate - batch_annotations.keys()
        else:
            logger.info('{} not using web'.format(self.name))

        return annotations
