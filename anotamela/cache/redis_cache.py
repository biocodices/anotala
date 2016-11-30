import logging

import redis

from anotamela.cache import Cache


logger = logging.getLogger(__name__)


class RedisCache(Cache):
    def __init__(self, host='localhost', port=6379, db=0):
        self.client = redis.StrictRedis(host=host, port=port, db=db)
        self.connection_data = self.client.connection_pool.connection_kwargs

        self.client.ping()  # Raises ConnectionError if server is not there
        logger.info('Connected to Redis @{host}:{port} db={db}'.format(
                **self.connection_data))

    def __repr__(self):
        rep = '{name}(host={host}, port={port}, db={db})'
        return rep.format({**self.connection_data,
                           'name': self.__class__.__name__})

    def _client_set(self, data_to_cache, namespace):
        data_to_cache = {self._id_to_key(id_, namespace): value
                         for id_, value in data_to_cache.items()}
        self.client.mset(data_to_cache)

    def _client_get(self, ids, namespace):
        keys = [self._id_to_key(id_, namespace) for id_ in ids]
        return dict(zip(ids, self.client.mget(keys)))

    @staticmethod
    def _id_to_key(id_, namespace):
        # Transform the ids to keys prefixing them with the given namespace
        return '{}:{}'.format(namespace, id_)

