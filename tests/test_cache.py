import logging

import pytest
import redis
from sqlalchemy.exc import OperationalError

from anotamela.cache import RedisCache, PostgresCache


logger = logging.getLogger()


TEST_INFO = {'key1': 'key1-val', 'key2': 'key2-val'}

def test_get(mock_cache):
    keys = 'key1 key2'.split()
    namespace = 'foo'
    mock_cache.storage[namespace].update(TEST_INFO)

    info_dict = mock_cache.get(keys, namespace)

    assert all([key in info_dict for key in keys])
    assert all([k in info_dict[k] for k in keys])

def test_set(mock_cache):
    mock_cache = mock_cache
    namespace = 'bar'

    mock_cache.set(TEST_INFO, namespace)
    cached_info = mock_cache.get(TEST_INFO.keys(), namespace)

    assert all(cached_info[k] == v for k, v in TEST_INFO.items())

def test_redis_cache():
    try:
        redis_cache = RedisCache()
        # FIXME: This is hacky. It will just test Redis if it finds it
        # in its default location (localhost:6379). Think of smth better.
    except redis.exceptions.ConnectionError:
        return

    namespace = get_test_namespace()

    redis_cache._client_set(TEST_INFO, namespace)
    cached_data = redis_cache.client.get('{}:key1'.format(namespace))

    assert cached_data == b'key1-val'

    cached_data = redis_cache.get(TEST_INFO.keys(), namespace)

    assert cached_data == TEST_INFO

    # Cleanup
    testing_keys = redis_cache.client.keys('{}*'.format(namespace))
    redis_cache.client.delete(*testing_keys)

def test_postgres_cache():
    try:
        postgres_cache = PostgresCache()
        # FIXME: This is hacky. It will just test Postgres if it finds it
        # in can connect with credentials found in the default path
        # ~/.postgres_credentials.yml. Ought to think of smth better.
    except OperationalError:
        return

    namespace = get_test_namespace()

    postgres_cache._client_set(TEST_INFO, namespace)
    cached_data = postgres_cache.get(TEST_INFO.keys(), namespace)

    assert cached_data == TEST_INFO

    # Cleanup
    postgres_cache._client_del(TEST_INFO.keys(), namespace)

def get_test_namespace():
    return 'testing_anotamela'

