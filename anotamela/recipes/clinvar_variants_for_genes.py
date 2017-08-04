from logging import getLogger

from Bio import Entrez
import pandas as pd

from anotamela.helpers import set_email_for_entrez
from anotamela.annotators import ClinvarVariationAnnotator


logger = getLogger(__name__)


COLUMN_ORDER = [
    'name',
    'chrom_g37',
    'start_g37',
    'stop_g37',
    'ref_g37',
    'alt_g37',
    'variation_type',
    'gene_symbol',
    'dbsnp_id',
    'associated_phenotypes',
    'clinical_significance',
    'genomic_change_g37',
    'coding_changes',
    'protein_changes',
    'consequences_functions',
    'frequencies',
    'omim_id',
    'uniprot_id',
    'variation_id',

    'genes',
    'clinical_summary',
    'clinical_assertions',
    'accession_g37',
    'length_g37',
    'allele_id',
    'alleles',
    'consequences',
    'genomic_change_g37_accession',
    'genomic_change_g37_name',
    'variant_type',
    'variation_name',
    'submitters',

    # GRCh38 fields
    'genomic_change_g38',
    'chrom_g38',
    'start_g38',
    'stop_g38',
    'ref_g38',
    'alt_g38',
    'accession_g38',
    'length_g38',
    'genomic_change_g38_accession',
    'genomic_change_g38_name',
]


def clinvar_variants_for_genes(gene_list, cache, as_dataframe=False):
    """Given a list of gene symbols (e.g. AIP, BRCA2), return a list of
    variants from a search in Clinvar with those gene symbols.

    Requires a *cache* argument for the ClinvarVariationAnnotator annotator.
    Options are: 'myqsl', 'postgres', 'dict', 'redis'.

    Returns a tuple of (variant_ids, annotations).

    If *as_dataframe* is set to True, it returns a pandas DataFrame with
    the annotations.
    """
    annotator = ClinvarVariationAnnotator(cache=cache)

    variant_ids = []
    variants = {}

    for gene_variant_ids in clinvar_variant_ids_for_genes(gene_list):
        variants_for_this_gene = annotator.annotate(gene_variant_ids)

        variant_ids += gene_variant_ids
        variants.update(variants_for_this_gene)

    if as_dataframe:
        df = pd.DataFrame(list(variants.values()))
        from pprint import pprint
        pprint(df.columns)
        print('-'*15)
        columns = [col for col in COLUMN_ORDER if col in df.columns]
        return df[columns]
    else:
        return (variant_ids, variants)

def clinvar_variant_ids_for_genes(gene_list):
    """Given a list of gene symbols, return a list of ClinVar variant IDs."""
    set_email_for_entrez()

    logger.info('Query ClinVar with {} genes: {}'.format(len(gene_list),
                                                         ', '.join(gene_list)))

    for gene in gene_list:
        query = '{}[Gene Name]'.format(gene)
        logger.info('ClinVar Query = "{}"'.format(query))
        handle = Entrez.esearch(db='clinvar', term=query, retmax=10000)
        # 10,000, big enough number to get every variant for a gene
        search_result = Entrez.read(handle)
        variant_ids = search_result['IdList']

        logger.info('Got {} ClinVar variant IDs for gene {}'
                    .format(len(variant_ids), gene))

        yield variant_ids

