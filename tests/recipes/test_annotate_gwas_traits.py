from anotamela import GwasCatalogAnnotator
from anotamela.recipes import annotate_gwas_traits


TEST_CACHE = 'mysql'


def test_annotate_gwas_traits(monkeypatch):
    def mock_annotate(self, ids, **kwargs):
        data = {
            'rs123': [{'trait': 'trait-123a'}, {'trait': 'trait-123b'}],
            'rs234': [{'trait': 'trait-234'}, {'trait': 'trait-234'}],
        }
        return {k: v for k, v in data.items() if k in ids}

    monkeypatch.setattr(GwasCatalogAnnotator, 'annotate', mock_annotate)

    result = annotate_gwas_traits(['rs123', 'rs234'], cache=TEST_CACHE)
    assert result == {'rs123': ['trait-123a', 'trait-123b'],
                      'rs234': ['trait-234']}

    result = annotate_gwas_traits('rs123', cache=TEST_CACHE)
    assert result == {'rs123': ['trait-123a', 'trait-123b']}

    assert annotate_gwas_traits('rs345', cache=TEST_CACHE) is None

