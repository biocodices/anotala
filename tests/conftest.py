pytest_plugins = ['helpers_namespace']

from os.path import join, dirname

import pytest
from sqlalchemy.exc import OperationalError

from anotamela.cache import create_cache
from anotamela.cache import AVAILABLE_CACHES


AVAILABLE_CACHES['mock_cache'] = AVAILABLE_CACHES['dict']


@pytest.fixture(scope='session')
def proxies():
    return {'http': 'socks5://localhost:9050'}


@pytest.helpers.register
def file(filename):
    return join(dirname(__file__), 'files', filename)


@pytest.helpers.register
def mock_annotate(ids, *args, **kwargs):
    fake_annotations = ['{}-ann'.format(id_) for id_ in ids]
    return dict(zip(ids, fake_annotations))


@pytest.helpers.register
def get_cache(cache_name):
    try:
        # Tests will write to the same database, but in a different table
        # (defined by the testing namespace of the TEST_PARAMS)
        return create_cache(cache_name)
    except OperationalError:
        return

