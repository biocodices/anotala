import pytest
from bs4 import BeautifulSoup

from anotamela.annotators import ClinvarVariationAnnotator


@pytest.fixture
def variations_xml():
    return """
<ClinvarResult-Set>
    <VariationReport VariationID="Var-1"
                     VariationName="Var-1-Name"
                     VariationType="Var-1-Type">

        <GeneList GeneCount="1">
            <Gene GeneID="Gene-1"
                  Symbol="Gene-1-Symbol"
                  Fullname="Gene 1 full name">
            </Gene>
        </GeneList>

    </VariationReport>
    <VariationReport VariationId="Var-2">
    </VariationReport>
<ClinvarResult-Set>
"""

@pytest.fixture
def variation_soup(variations_xml):
    soup = BeautifulSoup(variations_xml, 'lxml-xml')
    return soup.select_one('VariationReport')


def test_annotations_by_id(variations_xml):
    result = ClinvarVariationAnnotator._annotations_by_id(
        ids=['Var-1', 'Var-2'],
        xml=variations_xml,
    )
    result = dict(result)

    assert len(result) == 2
    assert result['Var-1'].startswith('<VariationReport VariationId="Var-1"')
    assert result['Var-2'].startswith('<VariationReport VariationId="Var-2"')


def test_extract_variation_id(variation_soup):
    result = ClinvarVariationAnnotator._extract_variation_id(variation_soup)
    assert result == 'Var-1'


def test_extract_variation_name(variation_soup):
    result = ClinvarVariationAnnotator._extract_variation_name(variation_soup)
    assert result == 'Var-1-Name'


def test_extract_variation_type(variation_soup):
    result = ClinvarVariationAnnotator._extract_variation_type(variation_soup)
    assert result == 'Var-1-Type'


def test_extract_genes(variation_soup):
    result = ClinvarVariationAnnotator._extract_genes(variation_soup)
    assert result == [
        {
            'id': 'Gene-1',
            'symbol': 'Gene-1-Symbol',
            'full_name': 'Gene 1 full name',
        }
    ]

