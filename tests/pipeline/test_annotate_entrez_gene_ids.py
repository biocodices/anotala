import pytest

from anotamela.annotators import ENTREZ_GENE_ANNOTATOR_CLASSES
from anotamela.pipeline import annotate_entrez_gene_ids


def test_annotate_entrez_gene_ids(monkeypatch):
    annotator_classes = ENTREZ_GENE_ANNOTATOR_CLASSES.values()

    for annotator_class in annotator_classes:
        monkeypatch.setattr(annotator_class, 'annotate',
                            pytest.helpers.mock_annotate)

    ids = [1, 2, 3]
    result = annotate_entrez_gene_ids(ids, cache='mock_cache', use_web=True,
                                      use_cache=True, proxies={'foo': 'bar'})

    # Check that the result dataframe is populated with all the annotators data
    for annotator_class in annotator_classes:
        assert annotator_class.SOURCE_NAME in result
