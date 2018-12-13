import pytest
from anotamela import GwasCatalogAnnotator


def test_url():
    assert GwasCatalogAnnotator._url('rs123') == \
        'https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/rs123'

    with pytest.raises(ValueError):
        GwasCatalogAnnotator._url('foobar')

def test_convert_keys_to_camelcase():
    GwasCatalogAnnotator.KEYS_TO_CAMELCASE = ['shouldChange']
    dic = {
        'shouldChange': 'foo',
        'shouldNotChange': 'bar',
    }
    result = GwasCatalogAnnotator._parse_annotation(dic)
    assert result['should_change'] == 'foo'
    assert result['shouldNotChange'] == 'bar'

#  def test_parse_annotation():
    #  raw = {
        #  'rsId': 'rs268',
        #  'functionalClass': 'missense_variant',
        #  'lastUpdateDate': '2018-01-01',
    #  }
    #  result = GwasCatalogAnnotator._parse_annotation(raw)
    #  print(result)
    #  assert result['rs_id'] == 'rs268'
    #  assert result['functional_class'] == 'missense'
    #  assert result['last_update_date'] == '2018-01-01'

# Data like:
# https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/rs7329174/associations?projection=associationBySnp
