from anotamela import SnpeffAnnotator


def test_snpeff_annotator_parse_hit():
    hit_single = {
        'allele': 'A',
        'snpeff': {'ann': {
            'effect': 'eff1&eff2',
            'hgvs_c': 'c.1A>T',
            'hgvs_p': 'p.Ala123*'
        }}
    }
    annotations = SnpeffAnnotator._parse_hit(hit_single)

    assert isinstance(annotations, list)
    assert len(annotations) == 1
    annotation = annotations[0]
    assert annotation['effects'] == ['eff1', 'eff2']
    assert 'genomic_allele' in annotation
    assert 'coding_allele' in annotation

    # Check prot change was parsed:
    assert annotation['hgvs_p'] == 'p.Ala123Ter'

    hit_multiple = {'allele': 'A',
                    'snpeff': {'ann': [{}, {}, {}]}}
    annotations = SnpeffAnnotator._parse_hit(hit_multiple)
    assert len(annotations) == 3

    hit_without_snpeff_data = {}
    assert SnpeffAnnotator._parse_hit(hit_without_snpeff_data) is None
