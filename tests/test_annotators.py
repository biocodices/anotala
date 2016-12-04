import pytest

from anotamela.annotators import (
        DbsnpWebAnnotator,
        DbsnpEntrezAnnotator,
        ClinvarRsAnnotator,
        HgvsAnnotator
    )


test_params = [
        (DbsnpWebAnnotator, {
            'ids_to_annotate': 'rs268 rs123',
            'keys_to_check': ('snp_id orgStr supported_assm is_clinical '
                              'assembly org')
        }),
        (DbsnpEntrezAnnotator, {
            'ids_to_annotate': 'rs268',
            'keys_to_check': ('hgvs alleles type links fxn '
                              'clinical_significance synonyms frequency')
        }),
        (ClinvarRsAnnotator, {
            'ids_to_annotate': 'rs268 rs199473059',
            'keys_to_check': ('alt gene rsid rcv type cytogenic hgvs '
                              'variant_id hg19 ref chrom hg38 allele_id')
        }),
        (HgvsAnnotator, {
            'ids_to_annotate': 'rs268 rs199473059',
            'keys_to_check': ('clinvar_hgvs_g clinvar_hgvs_c '
                              'myvariant_hgvs_g '
                              'snpeff_hgvs_c snpeff_hgvs_p ')

        })
    ]

@pytest.mark.parametrize('annotator_class,params', test_params)
def test_generic_annotator(annotator_class, params):
    ids_to_annotate = params['ids_to_annotate'].split()
    annotator = annotator_class(cache='mock_cache')

    # Test annotation from web
    info_dict = annotator.annotate(ids_to_annotate, use_cache=False)

    for id_ in ids_to_annotate:
        assert info_dict[id_]
        for key in params['keys_to_check'].split():
            print(id_, key)
            assert info_dict[id_][key]

    # Test the info was correctly cached after the first query
    cached_data = annotator.annotate(ids_to_annotate, use_web=False)
    for id_ in ids_to_annotate:
        assert cached_data[id_]

