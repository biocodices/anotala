from anotamela import EnsemblAnnotator


def test_parse_annotation():
    annotation = {
        'foo': 'bar',
        'populations': [],
        'population_genotypes': [],
        'genotypes': [],
    }
    parsed = EnsemblAnnotator._parse_annotation(annotation.copy())

    assert parsed['foo'] == 'bar'

    full_info_keys = [
        'populations',
        'population_genotypes',
        'genotypes'
    ]
    for key in full_info_keys:
        assert key not in parsed

    EnsemblAnnotator.full_info = True
    parsed = EnsemblAnnotator._parse_annotation(annotation.copy())

    for key in full_info_keys:
        assert key in parsed

