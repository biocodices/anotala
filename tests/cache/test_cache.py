import logging
import json

import pytest

from anotamela.cache import create_cache


logger = logging.getLogger()


TEST_PARAMS = [
    # (test_data, namespace, as_json)
    ({'key1': 'val1', 'key2': 'val2'}, '_anotamela_test', False),
    ({'key1': {'key1.1': 'val1.1',
               'key1.2': 'val1.2'}}, '_anotamela_test_json', True)
]


def test_get_nonexistent_key():
    cache = (pytest.helpers.get_cache('mysql') or
             pytest.helpers.get_cache('psql'))
    if cache:
        # Should not explode:
        cache.get('non-existent-key', '_anotamela_test', as_json=True)
        cache.get('non-existent-key', '_anotamela_test_json', as_json=False)


@pytest.mark.parametrize('test_data,namespace,as_json', TEST_PARAMS)
def test_get(namespace, test_data, as_json):
    mock_cache = create_cache('mock_cache')

    # Manually set the test values, manually json-dump if necessary:
    values_to_set = test_data
    if as_json:
        values_to_set = {k: json.dumps(v) for k, v in test_data.items()}
    mock_cache.storage[namespace].update(values_to_set)

    # Test the get method, which should automatically json-load if needed:
    test_keys = list(test_data.keys())
    cached_data = mock_cache.get(test_keys, namespace, as_json=as_json)
    assert all([cached_data[k] == test_data[k] for k in test_keys])


@pytest.mark.parametrize('test_data,namespace,as_json', TEST_PARAMS)
def test_set(namespace, test_data, as_json):
    mock_cache = create_cache('mock_cache')

    # Test the set method
    mock_cache.set(test_data, namespace, save_as_json=as_json)

    # Manually get the test values to compare them, json-load if necessary
    cached_data = mock_cache.storage[namespace]
    if as_json:
        cached_data = {k: json.loads(v) for k, v in cached_data.items()}
    assert all([cached_data[k] == test_data[k] for k in test_data])

