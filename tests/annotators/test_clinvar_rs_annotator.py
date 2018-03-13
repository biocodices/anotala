from anotamela import ClinvarRsAnnotator


def test_parse_annotation(monkeypatch):
    mock_parse_hit = lambda _: ['rcv-1', 'rcv-2']
    monkeypatch.setattr(ClinvarRsAnnotator, '_parse_hit', mock_parse_hit)

    # The mocked _parse_hit will transform each annotation {} into a list
    hits = [{'_id': 'A'}, {'_id': 'B'}, {'_id': 'C'}]

    # Make sure the lists are flattened, so we get 6 elements, not 3:
    assert len(ClinvarRsAnnotator._parse_annotation(hits)) == 6


def test_parse_hit_single_rcv():
    single_rcv_hit = {
        'allele': 'G',
        'clinvar': {
            'root_level_key': 'root_level_value',
            'variant_id': 123,
            'rcv': {
                'accession': 'RCV1',
                'clinical_significance': 'Pathogenic, risk factor',
                'preferred_name': 'NM_000237.2(LPL):c.953A>G (p.Asn318=)',
                'conditions': [{'identifiers': {}}],
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

    assert rcv_annotation['url']
    assert rcv_annotation['variant_url']
    assert rcv_annotation['genomic_allele']
    assert rcv_annotation['coding_allele']

    # Test prot change is parsed, with the '=' replace with an aminoacid
    assert rcv_annotation['prot_change'] == 'p.Asn318Asn'

    # Test identifiers are parsed
    condition = rcv_annotation['conditions'][0]
    assert 'identifiers' not in condition

    # Test clinsig is parsed
    assert 'clinical_significance' not in rcv_annotation
    assert rcv_annotation['clinical_significances'] == ['Pathogenic', 'risk factor']

    # Test hit info is copied to each RCV entry
    assert rcv_annotation['root_level_key'] == 'root_level_value'


def test_parse_preferred_name():
    func = ClinvarRsAnnotator._parse_preferred_name
    name = 'NM_000237.2(LPL):c.953A>G (p.Asn318Ser)'
    parsed = func(name)

    assert parsed['transcript'] == 'NM_000237.2'
    assert parsed['gene'] == 'LPL'
    assert parsed['cds_change'] == 'c.953A>G'
    assert parsed['prot_change'] == 'p.Asn318Ser'

    name = 'NM_000116.4(TAZ):c.646+14C>T'
    parsed = func(name)

    assert parsed['transcript'] == 'NM_000116.4'
    assert parsed['gene'] == 'TAZ'
    assert parsed['cds_change'] == 'c.646+14C>T'

    name = 'TTN:c.105180G>C (p.Glu35060Asp)'
    parsed = func(name)

    assert parsed['gene'] == 'TTN'
    assert parsed['cds_change'] == 'c.105180G>C'
    assert parsed['prot_change'] == 'p.Glu35060Asp'

    name = 'nothing to match'
    parsed = func(name)

    assert parsed == {}


def test_parse_hit_multiple_rcvs():
    multiple_rcv_hits = {
        'allele': 'A',
        'clinvar': {
            'variant_id': 1,
            'rcv': [{'conditions': {},  # dict instead of list!
                     'preferred_name': '',
                     'accession': '',
                     'clinical_significance': 'sig1'},

                    {'conditions': [],
                     'preferred_name': '',
                     'accession': '',
                     'clinical_significance': 'sig2'},

                    {'conditions': [],
                     'preferred_name': '',
                     'accession': '',
                     'clinical_significance': 'sig3, sig4'}]
        }
    }

    rcv_annotations = ClinvarRsAnnotator._parse_hit(multiple_rcv_hits)

    # One annotation per RCV entry:
    assert len(rcv_annotations) == 3

    # Listify the conditions
    assert all(isinstance(ann['conditions'], list) for ann in rcv_annotations)

    # Split significances
    assert [ann['clinical_significances'] for ann in rcv_annotations] == \
           [['sig1'], ['sig2'], ['sig3', 'sig4']]

    assert all('clinical_significance' not in ann for ann in rcv_annotations)


def test_parse_condition_identifiers():
    func = ClinvarRsAnnotator._parse_condition_identifiers
    assert func({}) == {}

    condition = {'identifiers': {'medgen': 'MG1', 'omim': 'MIM1'}}
    assert func(condition) == {'medgen_id': 'MG1', 'omim_id': 'MIM1'}

