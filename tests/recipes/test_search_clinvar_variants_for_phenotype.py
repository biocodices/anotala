from anotala.recipes.search_clinvar_variants_for_phenotypes import (
    _build_query_from_terms,
    _build_filters_from_clinsigs,
    _search_clinvar,
)


def test_build_query_from_terms():
    terms = ['Pheno-1', 'Pheno with spaces']
    query = _build_query_from_terms(terms)

    assert query == ('("Pheno-1"[Disease/Phenotype] OR '\
                     '"Pheno with spaces"[Disease/Phenotype])')


def test_build_filters_from_clinsigs():
    conditions = ['pathogenic', 'likely pathogenic']
    filters = _build_filters_from_clinsigs(conditions)

    assert filters == ('("clinsig pathogenic"[Properties] OR ' \
                       '"clinsig likely pathogenic"[Properties])')


def test_search_clinvar():
    query = _search_clinvar(pheno_terms=['Pheno-1'], clinsigs=['benign'],
                           dry_run=True)
    assert query == ('("Pheno-1"[Disease/Phenotype]) AND ' \
                     '("clinsig benign"[Properties])')

