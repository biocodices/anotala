import pytest
from bs4 import BeautifulSoup

from anotamela.annotators import ClinvarVariationAnnotator


@pytest.fixture
def variations_xml():
    # See https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&id=1550&rettype=variation
    # For an example of this kind of XML structure.
    return """
<ClinvarResult-Set>
    <VariationReport VariationID="Var-1"
                     VariationName="Var-1-Name"
                     VariationType="Var-1-Type">

        <GeneList GeneCount="1">
            <Gene GeneID="Gene-1"
                  Symbol="Gene-1-Symbol"
                  FullName="Gene 1 full name"
                  strand="-"
                  HGNCID="HGNC:123">

                <OMIM>234</OMIM>
            </Gene>
        </GeneList>

        <ClinicalAssertionList>
            <GermlineList>
                <Germline SubmitterName="Submitter-1"
                          DateLastSubmitted="2000-01-02">
                    <PhenotypeList>
                        <Phenotype Name="Pheno-1">
                            <XRefList>
                                <XRef ID="MIM-1" DB="OMIM" />
                            </XRefList>
                        </Phenotype>
                    </PhenotypeList>
                    <ClinicalSignificance>
                        <Description>ClinSig-1</Description>
                        <Method>Method-1</Method>
                    </ClinicalSignificance>
                </Germline>
            </GermlineList>
        </ClinicalAssertionList>

    </VariationReport>
    <VariationReport VariationID="Var-2">
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
    assert result['Var-1'].startswith('<VariationReport VariationID="Var-1"')
    assert result['Var-2'].startswith('<VariationReport VariationID="Var-2"')


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
            'entrez_id': 'Gene-1',
            'symbol': 'Gene-1-Symbol',
            'full_name': 'Gene 1 full name',
            'hgnc_id': '123',
            'strand': '-',
            'omim_id': '234',
        }
    ]

def test_extract_clinical_assertions(variation_soup):
    results = ClinvarVariationAnnotator._extract_clinical_assertions(variation_soup)
    assert len(results) == 1
    result = results[0]

    assert result == {
        'clinical_significance': 'ClinSig-1',
        'date_last_submitted': '2000-01-02',
        'method': 'Method-1',
        'phenotypes': [{'name': 'Pheno-1', 'omim_id': 'MIM-1'}],
        'submitter_name': 'Submitter-1',
    }

