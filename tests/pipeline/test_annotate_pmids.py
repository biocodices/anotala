import pytest

from anotamela.annotators import PubmedAnnotator
from anotamela.pipeline import annotate_pmids


def test_annotate_pmids(monkeypatch):
    monkeypatch.setattr(PubmedAnnotator, 'annotate',
                        pytest.helpers.mock_annotate)

    ids = ['PM1', 'PM1', 'PM1', 'PM2', 'PM3']
    result = annotate_pmids(ids, cache='mock_cache', use_web=True,
                            use_cache=True, proxies={'foo': 'bar'})

    assert all([id_ in result for id_ in ids])
    assert result['PM1'] == 'PM1-ann'

