from logging import getLogger

from Bio import Entrez

from anotamela.helpers import set_email_for_entrez
from anotamela.annotators import ClinvarVariationAnnotator


logger = getLogger(__name__)


def clinvar_variants_for_genes(gene_list, cache):
    """Given a list of gene symbols (e.g. AIP, BRCA2), return a list of
    variants from a search in Clinvar with those gene symbols.

    Requires a *cache* argument for the ClinvarVariationAnnotator annotator.
    """
    set_email_for_entrez()

    logger.info('Query ClinVar with genes: {}'.format(', '.join(gene_list)))

    query = ' OR '.join('{}[Gene Name]'.format(gene) for gene in gene_list)
    handle = Entrez.esearch(db='clinvar', term=query, retmax=10000) # Everything
    search_result = Entrez.read(handle)
    variant_ids = search_result['IdList']

    logger.info('Got {} ClinVar variant IDs'.format(len(variant_ids)))

    clinvar_variation = ClinvarVariationAnnotator(cache=cache)
    annotations = clinvar_variation.annotate(variant_ids)

    return annotations

