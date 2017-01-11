import pytest

from anotamela.cache import AVAILABLE_CACHES


AVAILABLE_CACHES['mock_cache'] = AVAILABLE_CACHES['_dict']


@pytest.fixture(scope='session')
def proxies():
    return {}

