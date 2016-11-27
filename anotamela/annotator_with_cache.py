from os import environ
import sys
import json
import redis
import time
import logging
from concurrent.futures import ThreadPoolExecutor

from tqdm import tqdm

from anotamela.helpers import grouped


logger = logging.getLogger(__name__)
coloredlogs.DEFAULT_LOG_FORMAT = '[@%(hostname)s %(asctime)s] %(message)s'
coloredlogs.install(level=environ.get('LOGLEVEL') or 'INFO')


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
    def __init__(self, redis_host='localhost', redis_port=6379, redis_db=0):
        self.name = self.__class__.__name__
        self.set_redis_client(redis_host, redis_port, redis_db)

    def set_redis_client(self, host, port=6379, db=0):
        # Set the Redis client as a class attribute so every instance shares
        # the same Redis connection:
        klass = self.__class__
        klass._redis_client = redis.StrictRedis(host=host, port=port, db=db)
        conn = klass._redis_client.connection_pool.connection_kwargs
        logger.info('{0} connected to Redis cache: {1}:{2} db={3}'.format(
                    self.name, conn['host'], conn['port'], conn['db']))

    def annotate(self, ids, parallel=10, sleep_time=10, use_cache=True,
                 use_web=True):
        """Annotate a list of variant identifiers. Annotations are cached."""
        if isinstance(ids, str) or isinstance(ids, int):
            ids = [ids]
        ids = set(ids)

        info_dict = {}

        if use_cache:
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

        self.parse_info_dict(info_dict)
        return info_dict

    def _batch_query(self, ids, parallel, sleep_time):
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
        return self._cache_get(ids, verbose=False)

    def _query_and_set_cache(self, id_):
        """Make a query for a single id, cache the response, and return it."""
        response = self._query(id_)
        if response:
            self._cache_set({id_: response})
        return response

    def _cache_set(self, info_dict):
        """
        Set the cache for a list of IDs. Expects a dict with the form:
        {
            id1: raw response of the service for ID 1,
            id2: raw response of the service for ID 2,
            ...
        }
        It will json-dump the dicts.
        """
        # Remove empty responses
        data_to_cache = {}
        for id_, value in info_dict.items():
            if not value:
                continue
            if isinstance(value, dict) or isinstance(value, list):
                value = json.dumps(value)
            data_to_cache[self._key(id_)] = value

        self._redis_client.mset(data_to_cache)

    def _cache_get(self, ids, verbose=True):
        """
        Get a list of IDs data from cache. Returns a dict with the form:
        {
            identifier1: <dict with info about the id1> or data_string,
            identifier2: <dict with info about id2> or data_string,
            ...
        }
        """
        keys = [self._key(id_) for id_ in ids]
        info_dict = {k: v for k, v in zip(ids, self._redis_client.mget(keys)) if v}

        for k, v in info_dict.items():
            v = self._decode_and_try_deserialize(v)
            info_dict[k] = v

        if verbose and len(ids) > 1:
            msg = '📂 %s found %s/%s in cache'
            logger.info(msg % (self.name, len(info_dict), len(ids)))

        return info_dict

    def _query(self):
        raise NotImplementedError()

    def _key(self):
        raise NotImplementedError()

    @staticmethod
    def _decode_and_try_deserialize(data):
        try:
            return json.loads(data.decode('utf8'))
        except ValueError:  # Not valid JSON, leave response as is
            return data.decode('utf8')

