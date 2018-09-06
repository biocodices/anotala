from anotamela.pipeline import generate_position_tags


def test_generate_position_tags():
    variant1 = {'dbsnp_web': {'rs_id': 'rs1',
                              'GRCh37.p13_chrom': '1',
                              'GRCh37.p13_start': 100,
                              'GRCh37.p13_stop': 100,
                              'GRCh38.p7_chrom': '1',
                              'GRCh38.p7_start': 110,
                              'GRCh38.p7_stop': 110}}
    result = generate_position_tags(variant1, assembly='GRCh37.p13')
    assert result == '1:100'

    variant2 = {'dbsnp_web': {'rs_id': 'rs2',
                              'GRCh37.p13_chrom': '2',
                              'GRCh37.p13_start': 200,
                              'GRCh37.p13_stop': 200,
                              'GRCh38.p7_chrom': '2',
                              'GRCh38.p7_start': 220,
                              'GRCh38.p7_stop': 220}}
    result = generate_position_tags(variant2, assembly='GRCh38.p7')
    assert result == '2:220'
