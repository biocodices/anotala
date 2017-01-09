from anotamela import DbnsfpAnnotator


def test_dbnsfp_annotator():
    info_1 = {'genename': 'Gene-1'}
    info_2 = {'genename': 'Gene-2'}

    raw = [
        {
            # This variant should be ommited, it doesn't have 'dbnsfp' data
            '_id': 'chr1:123A>G',
        },
        {
            '_id': 'chr1:123A>T',
            'dbnsfp': info_1,
        },
        {
            '_id': 'chr1:123A>C',
            'dbnsfp': info_2,
        }
    ]

    parsed = DbnsfpAnnotator._parse_annotation(raw)

    assert len(parsed) == 2
    assert all(info in parsed for info in [info_1, info_2])

