from collections import defaultdict

import pytest

from anotamela.cache import Cache
from anotamela.annotators.base_classes import AnnotatorWithCache


class MockCache(Cache):
    def __init__(self):
        self.name = 'Mock Cache'
        self.storage = defaultdict(dict)

    def _client_set(self, info_dict, namespace):
        table = self.storage[namespace]
        table.update(info_dict)

    def _client_get(self, keys, namespace):
        table = self.storage[namespace]
        return {k: table[k] for k in keys if k in table}

AnnotatorWithCache.AVAILABLE_CACHES['mock_cache'] = MockCache

@pytest.fixture(scope='module')
def mock_cache():
    return MockCache()

