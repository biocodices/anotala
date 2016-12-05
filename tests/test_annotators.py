import pytest

from anotamela.annotators import (
        DbsnpWebAnnotator,
        DbsnpEntrezAnnotator,
        DbsnpMyvariantAnnotator,
        ClinvarRsAnnotator,
        HgvsAnnotator,
        SnpeffAnnotator,
        MafAnnotator
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
        (DbsnpMyvariantAnnotator, {
            'ids_to_annotate': 'rs268 rs123',
            'keys_to_check': ('allele_origin alleles.allele alleles.freq '
                              'alt chrom class dbsnp_build flags gene.geneid '
                              'gene.symbol gmaf hg19.end hg19.start ref rsid '
                              'validated var_subtype vartype')
        }),
        (SnpeffAnnotator, {
            'ids_to_annotate': 'rs268 rs199473059',
            'keys_to_check': ('ann.feature_id ann.feature_type ann.gene_id '
                              'ann.genename ann.putative_impact ann.hgvs_c '
                              'ann.transcript_biotype')
        }),
        (ClinvarRsAnnotator, {
            'ids_to_annotate': 'rs268',
            'keys_to_check': ('accession allele_id alt cds_change chrom '
                              'clinical_significance conditions.medgen_id '
                              'conditions.name conditions.omim_id cytogenic '
                              'gene gene_id gene_symbol hg19_end hg19_start '
                              'hg38_end hg38_start hgvs_coding hgvs_genomic '
                              'last_evaluated number_submitters omim '
                              'origin preferred_name prot_change ref '
                              'review_status rsid transcript type url '
                              'variant_id')
        }),
        (HgvsAnnotator, {
            'ids_to_annotate': 'rs268 rs199473059',
            'keys_to_check': ('clinvar_hgvs_g clinvar_hgvs_c myvariant_hgvs_g '
                              'snpeff_hgvs_c snpeff_hgvs_p')
        }),
        (MafAnnotator, {
            'ids_to_annotate': 'rs268',
            'keys_to_check': ('1000gp3_sas_af 1000gp3_af exac_eas_af '
                              'exac_adj_af 1000gp3_gmaf exac_nfe_af exac_af '
                              'exac_amr_af 1000gp3_eas_af exac_afr_af '
                              'exac_fin_af 1000gp3_afr_af esp6500_aa_af '
                              'esp6500_ea_af exac_sas_af 1000gp3_amr_af '
                              '1000gp3_eur_af')
        })
    ]


@pytest.mark.parametrize('annotator_class,params', test_params)
def test_generic_annotator(annotator_class, params):
    ids_to_annotate = params['ids_to_annotate'].split()
    annotator = annotator_class(cache='mock_cache')

    # Test annotation from web
    info_dict = annotator.annotate(ids_to_annotate, use_cache=False)

    for id_ in ids_to_annotate:
        for key in params['keys_to_check'].split():
            check_dict_key(info_dict[id_], key)

    # Test the info was correctly cached after the first query
    cached_data = annotator.annotate(ids_to_annotate, use_web=False)
    for id_ in ids_to_annotate:
        for key in params['keys_to_check'].split():
            check_dict_key(cached_data[id_], key)


def check_dict_key(dictionary, key_to_check):
    # Some annotations are a list of dicts instead of a single dict
    # Handle those cases calling this function for each dict in the list:
    if isinstance(dictionary, list):
        for dict_item in dictionary:
            check_dict_key(dict_item, key_to_check)

    # If it's a dictionary
    else:
        if '.' in key_to_check:
            # If the key_to_check has a dot in it, split it and recall this
            # function going one level deeper, until there are no dots left
            top_key, sub_key = key_to_check.split('.', maxsplit=1)
            check_dict_key(dictionary[top_key], sub_key)
        else:
            assert dictionary[key_to_check] not in [None, [], {}]

