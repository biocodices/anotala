import pytest
from bs4 import BeautifulSoup

from anotamela.annotators import ClinvarRCVAnnotator


@pytest.fixture
def clinvar_set_xml():
    return """
<ClinVarSet ID="clinvarset-1">
    <ReferenceClinVarAssertion>
        <ClinVarAccession Acc="RCV-1" Type="RCV"/>
        <MeasureSet Type="Haplotype" ID="224889">
        </MeasureSet>
    </ReferenceClinVarAssertion>
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
    assert result['accession'] == 'RCV-1'
    assert result['entry_type'] == 'Haplotype'


def test_extract_accession(clinvar_set_soup):
    result = ClinvarRCVAnnotator._extract_accession(clinvar_set_soup)
    assert result == 'RCV-1'


def test_extract_entry_type(clinvar_set_soup):
    result = ClinvarRCVAnnotator._extract_entry_type(clinvar_set_soup)
    assert result == 'Haplotype'

