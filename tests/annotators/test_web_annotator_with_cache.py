import pytest

from anotamela.annotators.base_classes import WebAnnotatorWithCache


@pytest.fixture
def annotator():
    return WebAnnotatorWithCache('mock_cache')

def test_annotator_with_cache(annotator):
    inexistent_id = 'inexistent_id'
    assert annotator.annotate_one(inexistent_id, use_web=False) is None
    assert annotator.annotate(inexistent_id, use_web=False) == {}


def test_get_cached_ids(annotator):
    annotator.SOURCE_NAME = 'annotator-foo'
    annotator.cache.storage = {
        'annotator-foo': {
            'foo': 'foo!',
            'bar': 'bar!'
        }
    }
    assert annotator.get_cached_ids() == set(['foo', 'bar'])
