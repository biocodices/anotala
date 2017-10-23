from logging import getLogger

from Bio import Entrez

from anotamela.helpers import set_email_for_entrez
from anotamela.annotators import ClinvarVariationAnnotator


logger = getLogger(__name__)


def search_clinvar_variants_for_phenotypes(
        pheno_terms,
        clinsigs=['pathogenic', 'likely pathogenic'],
        cache='dict',
    ):
    """
    Queries ClinVar API for all the variants associated with any of the
    phenotypes passed as *pheno_terms*. Optionally, pass a list of *clinsigs*
    to filter the results. You may also pass a *cache* for
    ClinvarVariationAnnotator to use.

    Returns a list with those variants. Each variant is a dictionary with lots
    of data.

    Usage:

    ```
        from anotamela.recipes import search_clinvar_variants_for_phenotype

        phenos = ['Diabetes mellitus', 'Insulin']
        variants = search_clinvar_variants_for_phenotype(phenos)
    ```
    """
    variation_ids = _search_clinvar(pheno_terms, clinsigs)
    annotator = ClinvarVariationAnnotator(cache)
    annotations = annotator.annotate(variation_ids)
    return annotations


def _search_clinvar(pheno_terms, clinsigs, dry_run=False):
    """
    Take a list of phenotype/disease terms and a list of clinical significances
    and return a list of ClinVar Variation IDs that match the query.
    """
    pheno_query = _build_query_from_terms(pheno_terms)
    filters_query = _build_filters_from_clinsigs(clinsigs)
    query = ' AND '.join([pheno_query, filters_query])

    if dry_run:
        return query

    set_email_for_entrez()

    logger.info('Hit Entrez with query: "{}"'.format(query))
    handle = Entrez.esearch(db='clinvar', term=query, retmax=10000)
    response = Entrez.read(handle)
    ids = response['IdList']
    logger.info('Got {}/{} results'.format(len(ids), response['Count']))

    return ids


def _build_query_from_terms(pheno_terms):
    """
    Takes a list of phenotype terms like "Diabetes mellitus" and "Insulin"
    and returns an OR query.
    """
    pheno_queries = ['"{}"[Disease/Phenotype]'.format(term)
                     for term in pheno_terms]
    or_query = ' OR '.join(pheno_queries)
    return '({})'.format(or_query) # Surround in parens to avoid ambiguities


def _build_filters_from_clinsigs(clinsigs):
    """
    Takes a list of ClinVar clinical significances like "pathogenic" and
    "likely pathogenic" and returns an OR query of filters.
    """
    condition_queries = ['"clinsig {}"[Properties]'.format(clinsig)
                         for clinsig in clinsigs]
    or_query = ' OR '.join(condition_queries)
    return '({})'.format(or_query) # Surround in parens to avoid ambiguities

