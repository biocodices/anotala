import json
import pytest
import redis

from anotamela.cache import create_cache


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
        return create_cache('redis')
    except redis.exceptions.ConnectionError:
        return


def cleanup_redis_cache(cache, namespace):
    testing_keys = cache.client.keys('{}*'.format(namespace))
    cache.client.delete(*testing_keys)


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

