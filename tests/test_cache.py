import logging
import json

import pytest
import redis
from sqlalchemy.exc import OperationalError

from anotamela.cache import RedisCache, PostgresCache


logger = logging.getLogger()


TEST_PARAMS = [
    # (test_data, namespace, as_json)
    ({'key1': 'val1', 'key2': 'val2'}, '_anotamela_test', False),
    ({'key1': {'key1.1': 'val1.1',
               'key1.2': 'val1.2'}}, '_anotamela_test_json', True)
]


def get_redis_cache():
    # FIXME: This is hacky. It will just test Redis if it finds it
    # in the default location (localhost:6379). Think of smth better.
    try:
        return RedisCache()
    except redis.exceptions.ConnectionError:
        return


def cleanup_redis_cache(cache, namespace):
    testing_keys = cache.client.keys('{}*'.format(namespace))
    cache.client.delete(*testing_keys)


@pytest.mark.parametrize('test_data,namespace,as_json', TEST_PARAMS)
def test_get(mock_cache, namespace, test_data, as_json):
    # Manually set the test values, manually json-dump if necessary:
    values_to_set = test_data
    if as_json:
        values_to_set = {k: json.dumps(v) for k, v in test_data.items()}
    mock_cache.storage[namespace].update(values_to_set)

    # Test the get method, which should automatically json-load if needed:
    test_keys = list(test_data.keys())
    cached_data = mock_cache.get(test_keys, namespace, load_as_json=as_json)
    assert all([cached_data[k] == test_data[k] for k in test_keys])


@pytest.mark.parametrize('test_data,namespace,as_json', TEST_PARAMS)
def test_set(mock_cache, namespace, test_data, as_json):
    # Test the set method
    mock_cache.set(test_data, namespace, save_as_json=as_json)

    # Manually get the test values to compare them, json-load if necessary
    cached_data = mock_cache.storage[namespace]
    if as_json:
        cached_data = {k: json.loads(v) for k, v in cached_data.items()}
    assert all([cached_data[k] == test_data[k] for k in test_data])


@pytest.mark.parametrize('test_data,namespace,as_json', TEST_PARAMS)
def test_redis_set(namespace, test_data, as_json):
    redis_cache = get_redis_cache()
    if not redis_cache:
        return

    redis_cache._client_set(test_data, namespace, as_json)

    # Get that data manually and deserialize, to compare to original data
    test_keys = list(test_data.keys())
    first_key = test_keys[0]
    cached_value = redis_cache.client.get('{}:{}'.format(namespace, first_key))
    cached_value = cached_value.decode('utf-8')
    if as_json:
        cached_value = json.loads(cached_value)

    assert cached_value == test_data[first_key]
    cleanup_redis_cache(redis_cache, namespace)


@pytest.mark.parametrize('test_data,namespace,as_json', TEST_PARAMS)
def test_redis_get(namespace, test_data, as_json):
    redis_cache = get_redis_cache()
    if not redis_cache:
        return

    # Set data manually, manually json-dump if necessary
    test_keys = list(test_data.keys())
    values_to_set = {'{}:{}'.format(namespace, k): v
                     for k, v in test_data.items()}
    if as_json:
        values_to_set = {k: json.dumps(v) for k, v in values_to_set.items()}
    redis_cache.client.mset(values_to_set)

    # Test get retrieval and automatic deserialization when needed
    cached_data = redis_cache.get(test_keys, namespace, as_json)
    assert cached_data == test_data

    cleanup_redis_cache(redis_cache, namespace)


@pytest.mark.parametrize('test_data,namespace,as_json', TEST_PARAMS)
def test_postgres_cache(namespace, test_data, as_json):
    try:
        # Tests will write to the same database, but in a different table
        # (defined by the different namespace generated below)
        postgres_cache = PostgresCache('~/.postgres_credentials.yml')
    except OperationalError:
        return

    # Tests set and get
    postgres_cache._client_set(test_data, namespace, as_json)
    cached_data = postgres_cache.get(test_data.keys(), namespace, as_json)
    assert cached_data == test_data

    # Cleanup
    postgres_cache._client_del(test_data.keys(), namespace)

