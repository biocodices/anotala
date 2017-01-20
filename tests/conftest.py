pytest_plugins = ['helpers_namespace']

from os.path import join, dirname

import pytest

from anotamela.cache import AVAILABLE_CACHES


AVAILABLE_CACHES['mock_cache'] = AVAILABLE_CACHES['_dict']


@pytest.fixture(scope='session')
def proxies():
    return {'http': 'socks5://localhost:9050'}


@pytest.helpers.register
def file(filename):
    return join(dirname(__file__), 'files', filename)

