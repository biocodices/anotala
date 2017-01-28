from anotamela import DbsnpWebAnnotator


def test_parse_annotation():
    annotation = {
        'assembly': {
            'GRCh37.p13': [
                {'groupTerm': 'Primary_Assembly', 'snp2chrOrien': '1'},
                {'groupTerm': 'PATCHES', 'snp2chrOrien': '0'},
            ],

            'GRCh38.p7': [
                {'groupTerm': 'Primary_Assembly', 'snp2chrOrien': '0'},
                {'groupTerm': 'PATCHES', 'snp2chrOrien': '1'},
            ],
        }
    }

    parsed = DbsnpWebAnnotator._parse_annotation(annotation)
    # Explicitely test booleans, not truthiness:
    assert parsed['GRCh37.p13_reverse'] is True
    assert parsed['GRCh38.p7_reverse'] is False

