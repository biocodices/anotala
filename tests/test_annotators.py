import pytest

from anotamela.annotators import (
        DbsnpWebAnnotator,
        DbsnpEntrezAnnotator,
    )


keys_to_check_per_source = {
        'dbsnp_web': ('snp_id orgStr supported_assm is_clinical assembly org'),
        'dbsnp_entrez': ('hgvs alleles type links fxn clinical_significance '
                         'synonyms frequency'),
    }

@pytest.mark.parametrize('annotator_class,keys_to_check', [
        (DbsnpWebAnnotator, keys_to_check_per_source['dbsnp_web']),
        (DbsnpEntrezAnnotator, keys_to_check_per_source['dbsnp_entrez'])
    ])

def test_generic_annotator(annotator_class, keys_to_check):
    annotator = annotator_class(cache='mock_cache')
    ids_to_annotate = 'rs268 rs123'.split()
    info_dict =  annotator.annotate(ids_to_annotate)

    assert all(id_ in info_dict for id_ in ids_to_annotate)
    assert all(response for response in info_dict.values())  # No empty responses

    id_ = ids_to_annotate[0]
    for key in keys_to_check.split():
        assert info_dict[id_][key]

