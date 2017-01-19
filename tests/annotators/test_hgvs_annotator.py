from anotamela.annotators import HgvsAnnotator


def test_parse_hit():
    func = HgvsAnnotator._parse_hit

    hit = {'_id': 'chr1:g.123A>G'}
    parsed = func(hit)
    assert parsed['myvariant_hgvs_g'] == hit['_id']
    assert parsed['genomic_allele'] == 'G'

    hit['snpeff'] = {'ann': [{'hgvs_c': 'foo',
                              'hgvs_p': 'bar',
                              'feature_id': 'baz'}]}

    assert 'snpeff_hgvs_c' in func(hit)
    assert 'snpeff_hgvs_p' in func(hit)

    hit['clinvar'] = {'hgvs': {'genomic': ['spam'], 'coding': 'qux'}}
    assert 'clinvar_hgvs_g' in func(hit)
    assert 'clinvar_hgvs_c' in func(hit)


def test_parse_hgvs_from_clinvar():
    func = HgvsAnnotator._parse_hgvs_from_clinvar

    hit = {}
    assert func(hit) == {}

    hit['clinvar'] = {'hgvs': {'genomic': 'foo'}}
    assert func(hit)['clinvar_hgvs_g'] == 'foo'

    hit['clinvar']['hgvs']['coding'] = 'bar'
    assert func(hit)['clinvar_hgvs_c'] == 'bar'


def test_parse_hgvs_from_snpeff():
    func = HgvsAnnotator._parse_hgvs_from_snpeff

    hit = {}
    assert func(hit) == {}

    hit['snpeff'] = {'ann': {'feature_id': 'foo', 'hgvs_c': 'bar'}}
    assert func(hit)['snpeff_hgvs_c'] == ['foo:bar']

    hit['snpeff']['ann'] = [
        {'feature_id': 'feat-1', 'hgvs_c': 'cod-1', 'hgvs_p': 'prot-1'},
        {'feature_id': 'feat-2', 'hgvs_c': 'cod-2', 'hgvs_p': 'prot-2'},
    ]
    assert func(hit) == {
        'snpeff_hgvs_c': ['feat-1:cod-1', 'feat-2:cod-2'],
        'snpeff_hgvs_p': ['feat-1:prot-1', 'feat-2:prot-2']
    }

