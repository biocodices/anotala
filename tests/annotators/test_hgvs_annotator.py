from anotamela.annotators import HgvsAnnotator


def test_parse_hit():
    func = HgvsAnnotator._parse_hit

    hit = {'_id': 'chr1:g.123A>G'}
    parsed = func(hit)
    assert parsed['myvariant_hgvs_g'] == hit['_id']
    assert parsed['genomic_allele'] == 'G'

    hit['snpeff'] = {'ann': [{'hgvs_c': 'c.123A>G',
                              'hgvs_p': 'p.Ala123Gly',
                              'feature_id': 'NM_1.0'}]}

    parsed = func(hit)
    assert parsed['snpeff_hgvs_c'] == ['NM_1.0:c.123A>G']
    assert parsed['snpeff_hgvs_p'] == ['NM_1.0:p.Ala123Gly']

    hit['clinvar'] = {'hgvs': {'genomic': ['NG_1.0:c.124A>G'],
                               'coding': 'NM_1.0:c.124A>G'}}
    parsed = func(hit)
    assert parsed['clinvar_hgvs_g'] == ['NG_1.0:c.124A>G']
    assert parsed['clinvar_hgvs_c'] == 'NM_1.0:c.124A>G'

