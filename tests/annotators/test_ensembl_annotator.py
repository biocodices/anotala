import pandas as pd
from anotamela import EnsemblAnnotator


### FIXME uncomment !!!!!!!!!!!!!!!!!!!!!!!!!
#  def test_parse_annotation():
    #  annotation = {
        #  'MAF': '0.5',
        #  'foo': 'bar',
        #  'populations': [],
        #  'population_genotypes': [],
        #  'genotypes': [],
    #  }
    #  parsed = EnsemblAnnotator._parse_annotation(annotation.copy())

    #  assert parsed['foo'] == 'bar'
    #  assert parsed['MAF'] == 0.5

    #  full_info_keys = [
        #  'populations',
        #  'population_genotypes',
        #  'genotypes'
    #  ]
    #  for key in full_info_keys:
        #  assert key not in parsed

    #  EnsemblAnnotator.full_info = True
    #  parsed = EnsemblAnnotator._parse_annotation(annotation.copy())

    #  for key in full_info_keys:
        #  assert key in parsed


def test_parse_1KG_sample_name():
    result = EnsemblAnnotator.parse_1KG_sample_name('foo')
    assert result == 'foo'

    result = EnsemblAnnotator.parse_1KG_sample_name('1000GENOMES:foo:bar')
    assert result == 'bar'


def test_parse_genotype_string():
    result = EnsemblAnnotator.parse_genotype_string('T|C')
    assert result == ('T', 'C')

    result = EnsemblAnnotator.parse_genotype_string('T/A')
    assert result == ('T', 'A')

    result = EnsemblAnnotator.parse_genotype_string('G')
    assert result == ('G', )


def test_read_1kg_populations():
    result = EnsemblAnnotator.read_1kg_populations()
    assert len(result) == 2504
    assert list(result.columns) == ['sample', 'pop', 'super_pop', 'gender']
    assert result.loc[0, 'sample'] == 'HG00096'


def test_genotypes_list_to_dataframe():
    genotypes = [
        {
            'gender': 'Female',
            'sample': '1000GENOMES:phase_3:HG00097',
            'genotype': 'C|C'
        },
        {
            'gender': 'Male',
            'sample': '1000GENOMES:phase_3:HG02595',
            'genotype': 'T|T'
        },
        {
            'sample': 'Foo:Sample-1',
            'genotype': 'C'
        },
    ]

    result = EnsemblAnnotator.genotypes_list_to_dataframe(
        genotypes,
        variant_id='rs123',
        ref_allele='C'
    )
    assert len(result) == 3

    expected_columns = [
        'sample',
        'population',
        'region',
        'variant_id',
        'genotype_alleles',
        'ref_allele',
        'alt_allele_dosage',
        'in_1kg',
        'gender',
        'genotype',
        'sample_full_name',
    ]
    assert list(result.columns) == expected_columns
    assert list(result['in_1kg']) == [True, True, False]
    assert list(result['alt_allele_dosage']) == [2, 0, 1]
    assert list(result['variant_id']) == ['rs123', 'rs123', 'rs123']
    assert list(result['population']) == ['GBR', 'GWD', pd.np.nan]
    assert list(result['region']) == ['EUR', 'AFR', pd.np.nan]
