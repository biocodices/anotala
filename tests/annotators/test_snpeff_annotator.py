from anotamela import SnpeffAnnotator


def test_snpeff_annotator_parse_hit():
    hit_single = {
        'snpeff': {
            'ann': {
                'effect': 'eff-1&eff-2',
            }
        }
    }

    annotations = SnpeffAnnotator._parse_hit(hit_single)

    assert isinstance(annotations, list)
    assert len(annotations) == 1
    assert 'effect' not in annotations[0]
    assert annotations[0]['effects'] == ['eff-1', 'eff-2']

    hit_multiple = {
        'snpeff': {
            'ann': [
                {'effect': 'eff-3'},
                {'effect': 'eff-3'},
                {'effect': 'eff-4'},
            ]
        }
    }

    annotations = SnpeffAnnotator._parse_hit(hit_multiple)

    assert len(annotations) == 3
    assert annotations[0]['effects'] == ['eff-3']
    assert annotations[1]['effects'] == ['eff-3']
    assert annotations[2]['effects'] == ['eff-4']

