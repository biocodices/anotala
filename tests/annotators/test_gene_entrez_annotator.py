from anotamela import GeneEntrezAnnotator


def test_parse_annotation():
    raw_annotation = {
        'Chromosome': '8',
        'Description': 'Gene-Name',
        'MapLocation': 'Map-Loc',
        'Mim': ['mim-id'],
        'Name': 'Name-1',
        'NomenclatureName': 'Nomen-1',
        'NomenclatureSymbol': 'Symbol-1',
        'Organism': {'CommonName': 'Org-Name',
                     'ScientificName': 'Org-Sci-Name',
                     'TaxID': 'Tax-ID'},
        'OtherAliases': 'Alias-1, Alias-2',
        'OtherDesignations': 'Other-Gene-Name',
        'Summary': 'Summary-1',
        'id': '4023',
    }

    result = GeneEntrezAnnotator._parse_annotation(raw_annotation)

    assert(result['chromosome'] == '8')
    assert(result['description'] == 'Gene-Name')
    assert(result['map_location'] == 'Map-Loc')
    assert(result['mim'] == ['mim-id'])
    assert(result['name'] == 'Name-1')
    assert(result['nomenclature_name'] == 'Nomen-1')
    assert(result['nomenclature_symbol'] == 'Symbol-1')
    assert(result['organism_common_name'] == 'Org-Name')
    assert(result['organism_scientific_name'] == 'Org-Sci-Name')
    assert(result['organism_tax_id'] == 'Tax-ID')
    assert(result['other_aliases'] == 'Alias-1, Alias-2')
    assert(result['other_designations'] == 'Other-Gene-Name')
    assert(result['summary'] == 'Summary-1')
    assert(result['entrez_id'] == 4023)
    assert(result['url'] == 'https://www.ncbi.nlm.nih.gov/gene/4023')
