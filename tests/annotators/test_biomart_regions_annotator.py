from anotala.annotators import BiomartRegionsAnnotator


def test_parse_annotation():
    raw_response = ('dbSNP\trs1\t1\t100\t100\t1\n' +
                    'dbSNP\trs2\t1\t200\t200\t1\n' +
                    'dbSNP\trs3\t1\t300\t300\t1\n')

    parsed = BiomartRegionsAnnotator._parse_annotation(raw_response)
    assert len(parsed) == 3
    assert parsed[0] == {'source': 'dbSNP',
                         'rsid': 'rs1',
                         'chrom_g37': '1',
                         'chrom_start_g37': 100,
                         'chrom_end_g37': 100,
                         'chrom_strand_g37': 1}

