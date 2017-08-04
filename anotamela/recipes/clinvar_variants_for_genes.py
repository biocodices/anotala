from logging import getLogger

from Bio import Entrez

from anotamela.helpers import set_email_for_entrez
from anotamela.annotators import ClinvarVariationAnnotator


logger = getLogger(__name__)


def clinvar_variants_for_genes(gene_list, cache):
    """Given a list of gene symbols (e.g. AIP, BRCA2), return a list of
    variants from a search in Clinvar with those gene symbols.

    Requires a *cache* argument for the ClinvarVariationAnnotator annotator.
    Options are: 'myqsl', 'postgres', 'dict', 'redis'.

    Returns a tuple of (variant_ids, annotations).
    """
    annotator = ClinvarVariationAnnotator(cache=cache)

    variant_ids = []
    variants = {}

    for gene_variant_ids in clinvar_variant_ids_for_genes(gene_list):
        variants_for_this_gene = annotator.annotate(gene_variant_ids)

        variant_ids += gene_variant_ids
        variants.update(variants_for_this_gene)

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

