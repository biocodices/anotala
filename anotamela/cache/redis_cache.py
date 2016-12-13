import logging
import json

import redis

from anotamela.cache import Cache


logger = logging.getLogger(__name__)


class RedisCache(Cache):
    def __init__(self, host='localhost', port=6379, db=0):
        self.client = redis.StrictRedis(host=host, port=port, db=db)
        self.connection_data = self.client.connection_pool.connection_kwargs

        self.client.ping()  # Raises ConnectionError if server is not there
        logger.info('Connected to Redis ({host}:{port}/db={db})'.format(
                    **self.connection_data))

    def __repr__(self):
        rep = '{name}(host={host}, port={port}, db={db})'
        return rep.format({**self.connection_data,
                           'name': self.__class__.__name__})

    def _client_get(self, ids, namespace, load_as_json):
        keys = [self._id_to_key(id_, namespace) for id_ in ids]
        annotations = [ann.decode('utf-8') for ann in self.client.mget(keys)]
        # Redis client returns the values in the same order as the queried keys
        info_dict = dict(zip(ids, annotations))
        if load_as_json:
            info_dict = self._jsonload_dict_values(info_dict)
        return info_dict

    def _client_set(self, data_to_cache, namespace, save_as_json):
        data_to_cache = {self._id_to_key(id_, namespace): value
                         for id_, value in data_to_cache.items()}
        if save_as_json:
            data_to_cache = self._jsondump_dict_values(data_to_cache)
        self.client.mset(data_to_cache)

    @staticmethod
    def _id_to_key(id_, namespace):
        # Transform the ids to keys prefixing them with the given namespace
        return '{}:{}'.format(namespace, id_)

