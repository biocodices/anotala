import pytest

from anotamela.cache import Cache
from anotamela.annotators import AnnotatorWithCache


class MockCache(Cache):
    def __init__(self):
        self.name = 'Mock Cache'
        self.storage = {}

    def _client_set(self, info_dict):
        self.storage.update(info_dict)

    def _client_get(self, keys):
        return {k: self.storage[k] for k in keys if k in self.storage}

AnnotatorWithCache.AVAILABLE_CACHES['mock_cache'] = MockCache

@pytest.fixture(scope='module')
def mock_cache():
    return MockCache()

