import pytest

from anotamela import GwasCatalogAnnotator


def test_gwas_catalog_annotator_parse_annotation():
    groups = []
    raw = {'grouped': {'resourcename': {'groups': groups}}}

    assert GwasCatalogAnnotator._parse_annotation(raw) is None

    association_data = [{
        'rsId': ['rs1; rs2'],
        'strongestAllele': ['rs1-A; rs2-T'],
        'context': ['intron_variant; intron_variant'],
        'chromLocation': ['chr1:1; chr1:2'],
        'pValueExponent': -10,
        'pubmedId': 'PM1',
        'riskFrequency': '0.02',
        'betaDirection': 'NR',
        'numberOfIndividuals': [100, 50],
        'entrezMappedGenes': ['GENE1; GENE2'],
        'entrezMappedGeneLinks': ['furthest|2|-10|1',
                                  'mid|3|5|1',
                                  'closest|1|0|1'],
        'reportedGeneLinks': ['GENE1|1|ENSG01'],
        'ancestryLinks': ['initial|NR|U.S.|European|100|NR'],
    }]
    study_data = []
    disease_trait_data = []

    groups.extend([
        {
            'groupValue': 'association',
            'doclist': {'docs': association_data}
        },
        {
            'groupValue': 'study',
            'doclist': {'docs': study_data}
        },
        {
            'groupValue': 'diseasetrait',
            'doclist': {'docs': disease_trait_data}
        },
    ])

    parsed_data = GwasCatalogAnnotator._parse_annotation(raw)
    original = association_data[0]
    parsed = parsed_data[0]

    # Test it removes NR values
    assert 'beta_direction' not in parsed

    # Test it parses floats
    assert isinstance(parsed['risk_allele_frequency_in_controls'], float)

    # Test it parses colon separated values and associates the values with rsids
    multiannotation_fields = [
        'strongest_alleles',
        'allele_impacts',
        'chrom_locations',
        'entrez_mapped_gene_symbols',
    ]
    for field in multiannotation_fields:
        assert len(parsed[field]) == 2
        assert sorted(parsed[field]) == ['rs1', 'rs2']

    # Test it transforms the keys from camelCase to snake_case
    assert parsed['p_value_exponent'] == original['pValueExponent']

    # Test it renames some keys
    assert parsed['pmid'] == original['pubmedId']

    # Test it adds new data after parsing
    expected_new_keys = [
        'urls',
        'pubmed_entries',
        'entrez_mapped_genes',
        'reported_genes',
        'genomic_alleles',
        'sample_info',
    ]
    for key in expected_new_keys:
        assert key and parsed[key]

    # Test it sorts the mapped genes by relative position
    assert parsed['entrez_mapped_genes'][0]['symbol'] == 'closest'
    assert parsed['entrez_mapped_genes'][-1]['symbol'] == 'furthest'

    # Test it parses sample size
    assert 'initial_study' in parsed['sample_size']


@pytest.mark.parametrize('rsid_allele,expected_tuple', [
    ('rs123-A', 'A'),
    ('rs123-?', None),
])
def test_infer_allele(rsid_allele, expected_tuple):
    assert GwasCatalogAnnotator._infer_allele(rsid_allele) == expected_tuple


def test_infer_allele_raise():
    with pytest.raises(ValueError):
        GwasCatalogAnnotator._infer_allele('nonmatching-rsid')


@pytest.mark.parametrize('association_entry,expected_pubmed_entries', [
    ({'author_year_id': ['John Doe|2017|PM1'],
      'pmid': 'PM1'},
     [{'pmid': 'PM1', 'author': 'John Doe', 'year': '2017'}]),

    ({'author_year_id': ['John Doe|2017|PM1', 'Mary P|2017|PM3'],
      'pmid': 'PM2'},
     [{'pmid': 'PM1', 'author': 'John Doe', 'year': '2017'},
      {'pmid': 'PM2'},
      {'pmid': 'PM3', 'author': 'Mary P', 'year': '2017'}])
])
def test_parse_pubmed_entries(association_entry, expected_pubmed_entries):
    result = GwasCatalogAnnotator._parse_pubmed_entries(association_entry)
    assert result == expected_pubmed_entries


def test_parse_pubmed_entries_raise():
    with pytest.raises(Exception):
        association = {'author_year_id': ['John Doe|2017|PM1',
                                          'Johnny D|2017|PM1']}
        GwasCatalogAnnotator._parse_pubmed_entries(association)


@pytest.mark.parametrize('given_list,expected_dict', [
    ([100, 50], {'initial_study': 100, 'replication_study': 50}),
    ([100], {'initial_study': 100})
])
def test_parse_sample_size(given_list, expected_dict):
    result = GwasCatalogAnnotator._parse_sample_size(given_list)
    assert result == expected_dict


@pytest.mark.parametrize('gene_link,expected_result', [
    ('GENE1|1|0|X', {'symbol': 'GENE1',
                     'entrez_gene_id': '1',
                     'relative_position': 0,
                     'chromosome': 'X'}),
    ('GENE2|2|-10|19', {'symbol': 'GENE2',
                        'entrez_gene_id': '2',
                        'relative_position': -10,
                        'chromosome': '19'})
])
def test_parse_entrez_mapped_gene_link(gene_link, expected_result):
    result = GwasCatalogAnnotator._parse_entrez_mapped_gene_link(gene_link)
    assert result == expected_result


def test_parse_reported_gene_link():
    gene_link = 'GENE1|1|ENSG01'
    result = GwasCatalogAnnotator._parse_reported_gene_link(gene_link)
    expected_result = {'symbol': 'GENE1',
                       'entrez_gene_id': '1',
                       'ensembl_id': 'ENSG01'}
    assert result == expected_result


@pytest.mark.parametrize('ancestry_link,expected_info', [
    ('initial|NR|NR|European|100|NA',
     {'ancestries': ['European'],
      'sample_size': 100,
      'raw': 'initial|NR|NR|European|100|NA',
      'study_type': 'initial'}),

    ('replication|NR|Germany,U.K.|European|100|Augsburg, Germany;',
     {'ancestries': ['European'],
      'raw': 'replication|NR|Germany,U.K.|European|100|Augsburg, Germany;',
      'cities': ['Augsburg, Germany'],
      'countries': ['Germany', 'U.K.'],
      'sample_size': 100,
      'study_type': 'replication'})

])
def test_parse_ancestry_link(ancestry_link, expected_info):
    parsed_info = GwasCatalogAnnotator._parse_ancestry_link(ancestry_link)
    assert parsed_info == expected_info


def test_group_sample_info():
    initial_study = {'ancestries': ['Euro'], 'study_type': 'initial'}

    replication_studies = [
        {'ancestries': ['Afro'], 'study_type': 'replication'},
        {'ancestries': ['Asian'], 'study_type': 'replication'},
    ]

    sample_info = [initial_study] + replication_studies
    grouped = GwasCatalogAnnotator._group_sample_info(sample_info)
    assert grouped['initial'] == [initial_study]
    assert grouped['replication'] == replication_studies


@pytest.mark.parametrize('values,expected_result', [
    (['intron_variant'], ['intron_variant']),
    (['intron_variant; intron_variant'], ['intron_variant', 'intron_variant']),
    (['rs123-A; rs234-T'], ['rs123-A', 'rs234-T']),
    (['intron_variant; 3_prime_UTR_variant'],
     ['intron_variant', '3_prime_UTR_variant'])
])
def test_parse_colon_separated_values(values, expected_result):
    result = GwasCatalogAnnotator._parse_colon_separated_values(values)
    assert result == expected_result


def test_parse_ci_range():
    f = GwasCatalogAnnotator._parse_ci_range

    assert f('[NR]') is None
    assert f('[1.10-1.25]') == {'lower_limit': 1.1, 'upper_limit': 1.25}

