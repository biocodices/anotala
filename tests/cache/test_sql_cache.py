import pytest


TEST_PARAMS = [
    # (test_data, namespace, as_json)
    ({'key1': 'val1', 'key2': 'val2'}, '_anotamela_test', False),
    ({'key1': {'key1.1': 'val1.1',
               'key1.2': 'val1.2'}}, '_anotamela_test_json', True)
]


@pytest.mark.parametrize('test_data,namespace,as_json', TEST_PARAMS)
def test_cache(namespace, test_data, as_json):
    try:
        cache = pytest.helpers.get_cache('postgres')
    except OSError:  # config file not found
        pass
    else:
        _test_cache_operations(cache, test_data, namespace, as_json)

    try:
        cache = pytest.helpers.get_cache('mysql')
    except OSError:  # config file not found
        pass
    else:
        _test_cache_operations(cache, test_data, namespace, as_json)


def _test_cache_operations(cache_instance, test_data, namespace, as_json):
    # Tests set and get
    cache_instance._client_set(test_data, namespace, as_json)
    cached_data = cache_instance.get(test_data.keys(), namespace, as_json)
    assert cached_data == test_data

    # Cleanup
    cache_instance._client_del(test_data.keys(), namespace)

