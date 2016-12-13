import logging
import json

import pytest
import redis
from sqlalchemy.exc import OperationalError

from anotamela.cache import RedisCache, PostgresCache


logger = logging.getLogger()


TEST_PARAMS = [
    ('test_foo', {'key1': 'val1', 'key2': 'val2'}, False),
    ('test_bar', {'key1': {'key1.1': 'val1.1', 'key1.2': 'val1.2'}}, True)
]


@pytest.mark.parametrize('namespace,test_data,as_json', TEST_PARAMS)
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


@pytest.mark.parametrize('namespace,test_data,as_json', TEST_PARAMS)
def test_set(mock_cache, namespace, test_data, as_json):
    # Test the set method
    mock_cache.set(test_data, namespace, save_as_json=as_json)

    # Manually get the test values to compare them, json-load if necessary
    cached_data = mock_cache.storage[namespace]
    if as_json:
        cached_data = {k: json.loads(v) for k, v in cached_data.items()}
    assert all([cached_data[k] == test_data[k] for k in test_data])


@pytest.mark.parametrize('namespace,test_data,as_json', TEST_PARAMS)
def test_redis_cache(namespace, test_data, as_json):
    try:
        redis_cache = RedisCache()
        # FIXME: This is hacky. It will just test Redis if it finds it
        # in the default location (localhost:6379). Think of smth better.
    except redis.exceptions.ConnectionError:
        return

    # Test set
    redis_cache._client_set(test_data, namespace, as_json)
    first_key = list(test_data.keys())[0]
    cached_data = redis_cache.client.get('{}:{}'.format(namespace, first_key))
    if as_json:
        cached_data = json.loads(cached_data)
    assert cached_data == test_data[first_key]

    # Test get
    values_to_set = test_data
    if as_json:
        values_to_set = {k: json.dumps(v) for k, v in test_data.items()}
    redis_cache.client.mset(values_to_set)
    cached_data = redis_cache.get(test_data.keys(), namespace, as_json)
    assert cached_data == test_data

    # Cleanup
    testing_keys = redis_cache.client.keys('{}*'.format(namespace))
    redis_cache.client.delete(*testing_keys)


@pytest.mark.parametrize('namespace,test_data,as_json', TEST_PARAMS)
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

