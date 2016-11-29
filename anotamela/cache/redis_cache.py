import logging

import redis

from anotamela.cache import Cache


logger = logging.getLogger(__name__)


class RedisCache(Cache):
    def __init__(self, host='localhost', port=6379, db=0):
        self.client = redis.StrictRedis(host=host, port=port, db=db)
        self.name = self.__class__.__name__

        _connection_data = {**self.client.connection_pool.connection_kwargs,
                            'name': self.name}
        logger.info('{name} connected to Redis cache: {host}:{port} '
                    'db={db}'.format(**_connection_data))

    def _client_set(self, info_dict):
        self.client.mset(data_to_cache)

    def _client_get(self, keys):
        return {key: value for key, value in zip(keys, self.client.mget(keys))}

