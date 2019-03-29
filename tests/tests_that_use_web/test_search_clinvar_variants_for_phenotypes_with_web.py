from anotala.recipes import search_clinvar_variants_for_phenotypes


def test_search_clinvar_variants_for_phenotypes():
    # Smoke test just to check it's not broken in a basic way
    annotations = search_clinvar_variants_for_phenotypes(
        pheno_terms=['Diabetes mellitus', 'insulin'],
        clinsigs=['pathogenic'],
    )

    assert len(annotations) > 0
