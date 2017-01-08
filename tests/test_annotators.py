import pytest
import pandas as pd

from anotamela.annotators import *


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
        (DbnsfpAnnotator, {
            'ids_to_annotate': 'rs268',
            'keys_to_check': ('cds_strand gerp++ mutationassessor '
                              'mutationtaster polyphen2 provean sift')
        }
        ),
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
                              'phenotypes.id phenotypes.inheritance '
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
                              'Mesh Ids')
        }),
        (GwasCatalogAnnotator, {
            'ids_to_annotate': 'rs7329174',
            'keys_to_check': ('entrezMappedGenes ancestralGroups resourcename '
                              'publicationDate region riskFrequency '
                              'strongestAllele rsId pubmedId '
                              'entrezMappedGeneLinks chromosomePosition '
                              'publication multiSnpHaplotype chromosomeName '
                              'shortform_autosuggest synonym_autosuggest '
                              'replicateSampleDescription publicationLink '
                              'positionLinks context reportedGene '
                              'catalogPublishDate numberOfIndividuals '
                              'countriesOfRecruitment chromLocation author '
                              'label_autosuggest_ws shortForm merged '
                              'orPerCopyNum reportedGeneLinks id title '
                              'ancestryLinks label_autosuggest author_s '
                              'label synonym_autosuggest_e synonym traitName '
                              'initialSampleDescription traitName_s efoLink '
                              'pValueExponent currentSnp pValueMantissa '
                              'label_autosuggest_e mappedUri mappedLabel '
                              'locusDescription synonym_autosuggest_ws '
                              'platform snpInteraction parent range '
                              'accessionId _version_ studyId')
        })
    ]


@pytest.mark.parametrize('annotator_class,params', test_params)
def test_annotator(proxies, annotator_class, params):
    ids_to_annotate = params['ids_to_annotate'].split()
    annotator = annotator_class(cache='mock_cache')
    annotator.PROXIES = proxies

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

    # Test the annotators accept a single ID instead of a list of IDs
    id_ = params['ids_to_annotate'].split()[0]
    result = annotator.annotate(id_, use_web=False)
    assert id_ in result


def test_omim_annotate_from_entrez_ids(proxies):
    entrez_id = '63976'
    mim_id = '605557'
    gene_symbol = 'PRDM16'

    annotator = OmimGeneAnnotator(cache='mock_cache')
    annotator.PROXIES = proxies
    annotations = annotator.annotate_from_entrez_ids([entrez_id])
    variants = annotations[entrez_id]

    assert len(variants) == 6

    for variant in variants:
        assert variant['gene_omim_id'] == mim_id
        assert variant['gene_entrez_id'] == entrez_id
        assert variant['gene_symbol'] == gene_symbol


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

