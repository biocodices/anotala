import pytest
import pandas as pd

from anotamela.annotators import *
from anotamela.annotators.base_classes import EntrezAnnotator


test_params = [
        (DbsnpWebAnnotator, {
            'ids_to_annotate': 'rs268 rs123',
            'keys_to_check': ('snp_id orgStr supported_assm is_clinical '
                              'assembly org')
        }),
        (DbsnpEntrezAnnotator, {
            'ids_to_annotate': 'rs268',
            'keys_to_check': ('hgvs alleles type links fxn synonyms frequency')
        }),
        (DbsnpMyvariantAnnotator, {
            'ids_to_annotate': 'rs268 rs123',
            'keys_to_check': ('allele_origin alleles.allele alleles.freq '
                              'alt chrom class dbsnp_build flags gene.geneid '
                              'gene.symbol gmaf hg19.end hg19.start ref rsid '
                              'validated var_subtype vartype')
        }),
        (DbnsfpAnnotator, {
            'ids_to_annotate': 'rs268',
            'keys_to_check': ('cds_strand gerp++ mutationassessor '
                              'mutationtaster polyphen2 provean sift')
        }
        ),
        (SnpeffAnnotator, {
            'ids_to_annotate': 'rs268 rs199473059',
            'keys_to_check': ('feature_id feature_type gene_id genename '
                              'genomic_allele coding_allele putative_impact hgvs_c '
                              'transcript_biotype')
        }),
        (ClinvarRsAnnotator, {
            'ids_to_annotate': 'rs268',
            'keys_to_check': ('accession allele_id alt cds_change chrom '
                              'clinical_significances conditions.medgen_id '
                              'conditions.name conditions.omim_id cytogenic '
                              'gene gene_id gene_symbol hg19_end hg19_start '
                              'hg38_end hg38_start hgvs_coding hgvs_genomic '
                              'last_evaluated number_submitters omim '
                              'origin preferred_name prot_change ref '
                              'review_status rsid transcript type url '
                              'variant_id genomic_allele coding_allele')
        }),
        (ClinvarRCVAnnotator, {
            'ids_to_annotate': 'RCV000210468 RCV000148846',
            'keys_to_check': 'accession entry_type'
        }),
        (HgvsAnnotator, {
            'ids_to_annotate': 'rs268 rs199473059',
            'keys_to_check': ('clinvar_hgvs_g clinvar_hgvs_c myvariant_hgvs_g '
                              'snpeff_hgvs_c snpeff_hgvs_p')
        }),
        (FrequenciesAnnotator, {
            'ids_to_annotate': 'rs268',
            'keys_to_check': (
                'A.dbSNP.General|G.CADD_1000g.African (AFR)|G.dbNSFP_1000gp3.European (EUR)|'
                'G.dbNSFP_ESP6500|G.dbNSFP_ExAC.Finnish (FIN)|G.dbNSFP_twinsUK.General|'
                'G.dbSNP.General'
            ),
            'split_char': '|'
        }),
        (OmimGeneAnnotator, {
            'ids_to_annotate': '605557',
            'keys_to_check': ('gene_omim_id gene_name gene_symbol '
                              'linked_mim_ids phenotypes review_paragraphs '
                              'variant_id rsids prot_changes '
                              'gene_entrez_id incidental_gene '
                              'pubmed_entries.pmid pubmed_entries.citation '
                              'pubmed_entries.short_mention pubmed_entries.url')
        }),
        (OmimVariantAnnotator, {
            'ids_to_annotate': '605557.0003',
            'keys_to_check': ('gene_omim_id gene_name gene_symbol gene_url url '
                              'linked_mim_ids phenotypes.url '
                              'phenotypes.id phenotypes.inheritances '
                              'phenotypes.name prot_changes '
                              'pubmed_entries.pmid pubmed_entries.citation '
                              'pubmed_entries.short_mention pubmed_entries.url '
                              'review_paragraphs rsids variant_id '
                              'phenotypes.incidental')
        }),
        (MygeneAnnotator, {
            'ids_to_annotate': '4023',
            'keys_to_check': ('name symbol entrezgene MIM HGNC '
                              'swissprot type_of_gene'),
        }),
        (GeneEntrezAnnotator, {
            'ids_to_annotate': '4023',
            'keys_to_check': ('chromosome description name nomenclature_name '
                              'nomenclature_symbol organism_common_name '
                              'organism_scientific_name organism_tax_id '
                              'other_designations summary'),
        }),
        (UniprotAnnotator, {
            'ids_to_annotate': 'P06858',
            'keys_to_check': ('desc gene_symbol new_aa old_aa pmids pos '
                              'pos_stop prot_change prot_id prot_url review '
                              'rsid variant_id variant_url')
        }),
        (PubmedAnnotator, {
            'ids_to_annotate': '23788249',
            'keys_to_check': ('AMA_Citation CitationData Abstract ArticleDate '
                              'Mesh Ids pmid url')
        }),
        (GwasCatalogAnnotator, {
            'ids_to_annotate': 'rs7329174',
            'keys_to_check': ('rsids strongest_alleles p_value_exponent '
                              'sample_info sample_size entrez_mapped_genes '
                              'reported_genes pubmed_entries')
        })
    ]


@pytest.mark.parametrize('annotator_class,params', test_params)
def test_annotator(proxies, annotator_class, params):
    ids_to_annotate = params['ids_to_annotate'].split()
    annotator = annotator_class(cache='mock_cache', proxies=proxies)

    # Test annotation from web
    info_dict = annotator.annotate(ids_to_annotate, use_cache=False)

    split_char = params.get('split_char') or ' '

    for id_ in ids_to_annotate:
        for key in params['keys_to_check'].split(split_char):
            check_dict_key(info_dict[id_], key)

    # Test the info was correctly cached after the first query

    cached_data = annotator.annotate(ids_to_annotate, use_web=False)
    for id_ in ids_to_annotate:
        for key in params['keys_to_check'].split(split_char):
            check_dict_key(cached_data[id_], key)

    # Test the annotators accept a single ID instead of a list of IDs
    id_ = params['ids_to_annotate'].split()[0]
    result = annotator.annotate(id_, use_web=False)
    assert id_ in result

    # Test annotate_one
    assert result[id_] == annotator.annotate_one(id_, use_web=False)

    # Test the annotators can return the raw data instead of the parsed info
    if hasattr(annotator, '_parse_annotation'):
        raw = annotator.annotate(id_, use_web=False, parse=False)
        assert raw[id_]
        assert result[id_] != raw[id_]


def check_dict_key(dictionary, key_to_check):
    # Some annotations are a list of dicts instead of a single dict
    # Handle those cases calling this function for the first dict in the list:
    if isinstance(dictionary, list):
        check_dict_key(dictionary[0], key_to_check)
    # Some annotations are a pandas DataFrame. Check the key in each row:
    elif isinstance(dictionary, pd.DataFrame):
        for _, row in dictionary.iterrows():
            check_dict_key(row, key_to_check)
    # If it's a dictionary
    else:
        if '.' in key_to_check:
            # If the key_to_check has a dot in it, split it and recall this
            # function going one level deeper, until there are no dots left
            top_key, sub_key = key_to_check.split('.', maxsplit=1)
            check_dict_key(dictionary[top_key], sub_key)
        else:
            assert dictionary[key_to_check] is not None
            assert dictionary[key_to_check] is not []
            assert dictionary[key_to_check] is not {}

