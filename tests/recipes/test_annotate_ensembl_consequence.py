from anotamela import EnsemblAnnotator
from anotamela.recipes import annotate_ensembl_consequence


def test_annotate_ensembl_consequence(monkeypatch):
    def mock_annotate(self, ids, **kwargs):
        data = {
            'rs123': {'most_severe_consequence': 'conseq-123'},
            'rs234': {'most_severe_consequence': 'conseq-234'},
        }
        return {k: v for k, v in data.items() if k in ids}

    monkeypatch.setattr(EnsemblAnnotator, 'annotate', mock_annotate)

    result = annotate_ensembl_consequence(['rs123', 'rs234'])
    assert result == {'rs123': 'conseq-123', 'rs234': 'conseq-234'}

    result = annotate_ensembl_consequence('rs123')
    assert result == {'rs123': 'conseq-123'}

