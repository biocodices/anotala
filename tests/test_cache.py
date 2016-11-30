import pytest


TEST_INFO = {'key1': 'key1-val', 'key2': 'key2-val'}

def test_get(mock_cache):
    keys = 'key1 key2'.split()
    namespace = 'foo'
    mock_cache.storage[namespace].update(TEST_INFO)

    info_dict = mock_cache.get(keys, namespace=namespace)

    assert all([key in info_dict for key in keys])
    assert all([k in info_dict[k] for k in keys])

def test_set(mock_cache):
    mock_cache = mock_cache
    namespace = 'bar'

    mock_cache.set(TEST_INFO, namespace=namespace)
    cached_info = mock_cache.get(TEST_INFO.keys(), namespace=namespace)

    assert all(cached_info[k] == v for k, v in TEST_INFO.items())

