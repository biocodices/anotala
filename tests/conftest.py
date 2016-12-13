from collections import defaultdict
import json

import pytest

from anotamela.cache import Cache
from anotamela.annotators.base_classes import AnnotatorWithCache


class MockCache(Cache):
    def __init__(self):
        self.name = 'Mock Cache'
        self.storage = defaultdict(dict)

    def _client_set(self, info_dict, namespace, save_as_json):
        table = self.storage[namespace]
        table.update(info_dict)

    def _client_get(self, keys, namespace, load_as_json):
        table = self.storage[namespace]
        info_dict = {k: table[k] for k in keys if k in table}
        if load_as_json:
            info_dict = {k: json.loads(v) for k, v in info_dict.items()}
        return info_dict

AnnotatorWithCache.AVAILABLE_CACHES['mock_cache'] = MockCache

@pytest.fixture(scope='module')
def mock_cache():
    return MockCache()

