from anotamela.pipeline import extract_entrez_genes


def test_extract_entrez_genes():
    dbsnp_entries = [
        {'gene': [{'geneid': 1, 'symbol': 'GENE-1'}]},
        {'gene': [{'geneid': 2, 'symbol': 'GENE-2'}]},
    ]

    ids = extract_entrez_genes(dbsnp_entries, field='geneid')
    assert set(ids) == {1, 2}

    symbols = extract_entrez_genes(dbsnp_entries, field='symbol')
    assert set(symbols) == {'GENE-1', 'GENE-2'}

