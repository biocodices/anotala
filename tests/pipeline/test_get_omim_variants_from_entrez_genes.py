from anotamela import OmimGeneAnnotator
from anotamela.pipeline import get_omim_variants_from_entrez_genes


def test_get_omim_variants_from_entrez_genes(monkeypatch):
    def mock_annotate_from_entrez_ids(self, entrez_ids, **kwargs):
        return {entrez_id: [{'foo': 'bar'}, {'baz': 'qux'}]
                for entrez_id in entrez_ids}

    monkeypatch.setattr(OmimGeneAnnotator, 'annotate_from_entrez_ids',
                        mock_annotate_from_entrez_ids)

    omim_variants = get_omim_variants_from_entrez_genes(['123', '234'],
                                                        cache='mock_cache')

    # Test the result is a flat list of entries:
    assert all(isinstance(entry, dict) for entry in omim_variants)

    # Test all the variants were included (2 variants * 2 annotated entrez ids)
    assert len(omim_variants) == 4

