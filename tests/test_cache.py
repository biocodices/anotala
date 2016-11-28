import pytest


from anotamela.cache import Cache


def test_get():
    generic_cache = Cache()

    def mock_client_get(keys):
        return ['{}-val'.format(key) for key in keys]

    generic_cache.name = 'Generic Cache'
    generic_cache._client_get = mock_client_get

    keys = 'key1 key2'.split()
    info_dict = generic_cache.get(keys)

    assert all(key in info_dict for key in keys)
    assert all(k in info_dict[k] for k in keys)

def test_set():
    generic_cache = Cache()

    cached_data = {}
    def mock_client_set(parsed_info_dict):
        nonlocal cached_data
        cached_data.update(parsed_info_dict)

    generic_cache.name = 'Generic Cache'
    generic_cache._client_set = mock_client_set

    info_dict = {'key1': 'key1-val', 'key2': 'key2-val'}
    generic_cache.set(info_dict)

    assert all(cached_data[k] == v for k, v in info_dict.items())

