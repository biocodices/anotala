from anotamela import ClinvarRsAnnotator


def test_clinvar_annotator_parse_annotation(monkeypatch):
    mock_parse_hit = lambda _: ['rcv-1', 'rcv-2']
    monkeypatch.setattr(ClinvarRsAnnotator, '_parse_hit', mock_parse_hit)

    # The mocked _parse_hit will transform each annotation {} into a list
    hits = [{}, {}, {}]

    # Make sure the lists are flattened, so we get 6 elements, not 3:
    assert len(ClinvarRsAnnotator._parse_annotation(hits)) == 6


def test_clinvar_annotator_parse_hit_single_rcv():
    single_rcv_hit = {
        'clinvar': {
            'root_level_key': 'root_level_value',
            'rcv': {
                'accession': 'RCV1',
                'clinical_significance': 'Pathogenic, risk factor',
                'preferred_name': 'NM_000237.2(LPL):c.953A>G (p.Asn318Ser)',
                'conditions': [{
                    'identifiers': {'omim': 'MIM-1', 'medgen': 'MG-1'}
                }]
            }
        }
    }

    rcv_annotations = ClinvarRsAnnotator._parse_hit(single_rcv_hit)

    # Test it returns a list of RCV entries
    assert isinstance(rcv_annotations, list)
    assert len(rcv_annotations) == 1

    rcv_annotation = rcv_annotations[0]

    # Test that it parses the preferred name
    assert rcv_annotation['cds_change'] == 'c.953A>G'

    assert 'url' in rcv_annotation

    # Test identifiers are parsed
    condition = rcv_annotation['conditions'][0]
    assert condition['omim_id'] == 'MIM-1'
    assert condition['medgen_id'] == 'MG-1'
    assert 'identifiers' not in condition

    # Test clinsig is parsed
    assert 'clinical_significance' not in rcv_annotation
    assert rcv_annotation['clinical_significances'] == ['Pathogenic', 'risk factor']

    # Test hit info is copied to each RCV entry
    assert rcv_annotation['root_level_key'] == 'root_level_value'


def test_clinvar_annotator_parse_preferred_name():
    name = 'NM_000237.2(LPL):c.953A>G (p.Asn318Ser)'
    parsed = ClinvarRsAnnotator._parse_preferred_name(name)

    assert parsed['transcript'] == 'NM_000237.2'
    assert parsed['gene'] == 'LPL'
    assert parsed['cds_change'] == 'c.953A>G'
    assert parsed['prot_change'] == 'p.Asn318Ser'

    name = 'NM_000116.4(TAZ):c.646+14C>T'
    parsed = ClinvarRsAnnotator._parse_preferred_name(name)

    assert parsed['transcript'] == 'NM_000116.4'
    assert parsed['gene'] == 'TAZ'
    assert parsed['cds_change'] == 'c.646+14C>T'

    name = 'TTN:c.105180G>C (p.Glu35060Asp)'
    parsed = ClinvarRsAnnotator._parse_preferred_name(name)

    assert parsed['gene'] == 'TTN'
    assert parsed['cds_change'] == 'c.105180G>C'
    assert parsed['prot_change'] == 'p.Glu35060Asp'

    name = 'nothing to match'
    parsed = ClinvarRsAnnotator._parse_preferred_name(name)

    assert parsed == {}


def test_clinvar_annotator_parse_hit_multiple_rcvs():
    multiple_rcv_hits = {
        'clinvar': {
            'rcv': [{'conditions': [],
                     'preferred_name': '',
                     'accession': '',
                     'clinical_significance': ''},

                    {'conditions': [],
                     'preferred_name': '',
                     'accession': '',
                     'clinical_significance': ''},

                    {'conditions': [],
                     'preferred_name': '',
                     'accession': '',
                     'clinical_significance': ''}]
        }
    }

    rcv_annotations = ClinvarRsAnnotator._parse_hit(multiple_rcv_hits)
    assert len(rcv_annotations) == 3

