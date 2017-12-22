import pytest
from bs4 import BeautifulSoup

from anotamela.annotators import ClinvarRCVAnnotator


@pytest.fixture
def clinvar_set_xml():
    return """
<ClinVarSet ID="clinvarset-1">
    <Title>
        Title-1
    </Title>
    <ReferenceClinVarAssertion>
        <ClinVarAccession Acc="RCV-1" Type="RCV"/>
        <MeasureSet Type="Haplotype" ID="224889">
        </MeasureSet>
    </ReferenceClinVarAssertion>
    <MeasureSet ID="123" Type="Type-1">
        <Measure ID="234" Type="Type-2">
            <XRef DB="dbSNP" ID="123" Type="rs" />
            <AttributeSet>
                <Attribute Accession="Acc-1" Change="Change-1" Type="Type-1">
                    Full-Name-1
                </Attribute>
            </AttributeSet>
            <AttributeSet>
                <Attribute Accession="Acc-2" Change="Change-2" Type="Type-2">
                    Full-Name-2
                </Attribute>
            </AttributeSet>
        </Measure>
    </MeasureSet>
</ClinVarSet>
"""

@pytest.fixture
def clinvar_set_soup(clinvar_set_xml):
    return BeautifulSoup(clinvar_set_xml, 'lxml-xml')


def test_annotations_by_id(clinvar_set_xml):
    multi_accession_xml = """<?xml version="1.0" encoding="UTF-8" ?>
<ClinVarResult-Set>
    {}
    <ClinVarSet ID="clinvarset-2">
        <ReferenceClinVarAssertion>
            <ClinVarAccession Acc="RCV-2" Type="RCV"/>
        </ReferenceClinVarAssertion>
    </ClinVarSet>
</ClinVarResult-Set>
""".format(clinvar_set_xml)

    result = ClinvarRCVAnnotator._annotations_by_id([], multi_accession_xml)
    result = dict(result)

    assert 'clinvarset-1' in result['RCV-1']
    assert 'clinvarset-2' not in result['RCV-1']
    assert 'clinvarset-2' in result['RCV-2']
    assert 'clinvarset-1' not in result['RCV-2']

    assert result['RCV-1'].startswith('<ClinVarSet')
    assert result['RCV-2'].startswith('<ClinVarSet')


def test_parse_annotation(clinvar_set_xml):
    result = ClinvarRCVAnnotator._parse_annotation(clinvar_set_xml)
    keys = [
        'accession',
        'entry_type',
        'title',
        'attributes',
        'dbsnp_id'
    ]
    for key in keys:
        assert key in result


def test_extract_title(clinvar_set_soup):
    result = ClinvarRCVAnnotator._extract_title(clinvar_set_soup)
    assert result == 'Title-1'


def test_extract_attributes(clinvar_set_soup):
    result = ClinvarRCVAnnotator._extract_attributes(clinvar_set_soup)
    assert result == [
        {
            'accession': 'Acc-1',
            'change': 'Change-1',
            'type': 'Type-1',
            'full_name': 'Full-Name-1'
        },
        {
            'accession': 'Acc-2',
            'change': 'Change-2',
            'type': 'Type-2',
            'full_name': 'Full-Name-2'
        }
    ]


def test_extract_accession(clinvar_set_soup):
    result = ClinvarRCVAnnotator._extract_accession(clinvar_set_soup)
    assert result == 'RCV-1'


def test_extract_dbsnp_id(clinvar_set_soup):
    result = ClinvarRCVAnnotator._extract_dbsnp_id(clinvar_set_soup)
    assert result == 'rs123'


def test_extract_entry_type(clinvar_set_soup):
    result = ClinvarRCVAnnotator._extract_entry_type(clinvar_set_soup)
    assert result == 'Haplotype'

