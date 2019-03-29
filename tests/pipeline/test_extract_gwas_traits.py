from anotala.pipeline import extract_gwas_traits


def test_extract_gwas_traits():
    gwas_annotations = [
        {'trait': 'trait-2'},
        {'trait': 'trait-1'},
        {'trait': 'trait-3'},
        {'trait': 'trait-1'},
        {}
    ]

    result = extract_gwas_traits(gwas_annotations)

    assert result == ['trait-1', 'trait-2', 'trait-3']

