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

        <Allele AlleleID="Allele-1">
            <Name>Allele-Name-1</Name>
            <VariantType></VariantType>
        </Allele>

    </VariationReport>
</ClinvarResult-Set>
"""

@pytest.fixture
def variation_soup(variations_xml):
    soup = BeautifulSoup(variations_xml, 'lxml-xml')
    return soup.select_one('VariationReport')

def make_soup(xml):
    return BeautifulSoup(xml, 'lxml-xml').findChildren()[0]


def test_annotations_by_id():
    variations_xml = """
        <ClinvarResult-Set>
            <VariationReport VariationID="Var-1"></VariationReport>
            <VariationReport VariationID="Var-2"></VariationReport>
        </ClinvarResult-Set>
    """
    result = ClinvarVariationAnnotator._annotations_by_id(
        ids=['Var-1', 'Var-2'],
        xml=variations_xml,
    )
    result = dict(result)

    assert len(result) == 2
    assert result['Var-1'].startswith('<VariationReport VariationID="Var-1"')
    assert result['Var-2'].startswith('<VariationReport VariationID="Var-2"')


def test_extract_variation_id():
    soup = make_soup('<VariationReport VariationID="Var-1"></VariationReport')
    result = ClinvarVariationAnnotator._extract_variation_id(soup)
    assert result == 'Var-1'


def test_extract_variation_name():
    soup = make_soup('<VariationReport VariationName="Var-1-Name"></VariationReport')
    result = ClinvarVariationAnnotator._extract_variation_name(soup)
    assert result == 'Var-1-Name'


def test_extract_variation_type():
    soup = make_soup('<VariationReport VariationType="Var-1-Type"></VariationReport')
    result = ClinvarVariationAnnotator._extract_variation_type(soup)
    assert result == 'Var-1-Type'


def test_extract_genes():
    soup = make_soup("""
        <VariationReport>
            <GeneList GeneCount="1">
                <Gene GeneID="Gene-1"
                        Symbol="Gene-1-Symbol"
                        FullName="Gene 1 full name"
                        strand="-"
                        HGNCID="HGNC:123">

                    <OMIM>234</OMIM>
                </Gene>
            </GeneList>
        </VariationReport>
    """)
    result = ClinvarVariationAnnotator._extract_genes(soup)
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


def test_extract_clinical_assertions():
    soup = make_soup("""
        <VariationReport>
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
    """)
    results = ClinvarVariationAnnotator._extract_clinical_assertions(soup)
    assert len(results) == 1
    result = results[0]

    assert result == {
        'clinical_significance': 'ClinSig-1',
        'date_last_submitted': '2000-01-02',
        'method': 'Method-1',
        'phenotypes': [{'name': 'Pheno-1', 'omim_id': 'MIM-1'}],
        'submitter_name': 'Submitter-1',
    }


def test_extract_allele_basic_info():
    soup = make_soup("""
        <Allele AlleleID="Allele-1">
            <Name>Allele-Name</Name>
            <VariantType>Var-Type</VariantType>
        </Allele>
    """)

    result = ClinvarVariationAnnotator._extract_allele_basic_info(soup)

    assert result['name'] == 'Allele-Name'
    assert result['allele_id'] == 'Allele-1'
    assert result['variant_type'] == 'Var-Type'


def test_extract_sequence_info_from_allele():
    soup = make_soup("""
        <Allele>
            <SequenceLocation Assembly="GRCh37"
                                Chr="X"
                                Accession="NC_1"
                                start="1000"
                                stop="1000"
                                variantLength="1"
                                referenceAllele="A"
                                alternateAllele="G"/>
            <SequenceLocation Assembly="GRCh38"
                                Chr="X"
                                Accession="NC_2"
                                start="2000"
                                stop="2000"
                                variantLength="1"
                                referenceAllele="A"
                                alternateAllele="G"/>
        </Allele>
    """)

    result = ClinvarVariationAnnotator._extract_sequence_info_from_allele(soup)

    assert result['start_g37'] == 1000
    assert result['stop_g37'] == 1000
    assert result['chrom_g37'] == 'X'
    assert result['length_g37'] == 1
    assert result['ref_g37'] == 'A'
    assert result['alt_g37'] == 'G'
    assert result['accession_g37'] == 'NC_1'

    assert result['start_g38'] == 2000
    assert result['stop_g38'] == 2000
    assert result['chrom_g38'] == 'X'
    assert result['length_g38'] == 1
    assert result['ref_g38'] == 'A'
    assert result['alt_g38'] == 'G'
    assert result['accession_g38'] == 'NC_2'


def test_extract_allele_hgvs():
    soup = make_soup("""
        <Allele>
            <HGVSlist>
                <HGVS Assembly="GRCh37"
                      Change="genomic-change-37"
                      AccessionVersion="Accession-1"></HGVS>
                <HGVS Assembly="GRCh38"
                      Change="genomic-change-38"
                      AccessionVersion="Accession-1"></HGVS>
                <HGVS Type="HGVS, coding, RefSeq"
                      Change="coding-change"
                      AccessionVersion="Accession-1"></HGVS>
                <HGVS Type="HGVS, protein, RefSeq"
                      Change="protein-change"
                      AccessionVersion="Accession-1"></HGVS>
            </HGVSlist>
        </Allele>
    """)

    result = ClinvarVariationAnnotator._extract_allele_hgvs(soup)

    assert result['genomic_change_g37'] == 'genomic-change-37'
    assert result['genomic_change_g37_accession'] == 'Accession-1'

    assert result['genomic_change_g38'] == 'genomic-change-38'
    assert result['genomic_change_g38_accession'] == 'Accession-1'

    assert result['coding_change'] == 'coding-change'
    assert result['coding_change_accession'] == 'Accession-1'

    assert result['protein_change'] == 'protein-change'
    assert result['protein_change_accession'] == 'Accession-1'


def test_extract_xrefs():
    soup = make_soup("""
        <Allele>
            <XRefList>
                <XRef DB="UniProtKB" ID="Uniprot-1" />
                <XRef DB="OMIM" ID="Omim-1" />
                <XRef DB="dbSNP" ID="123" Type="rs" />
            </XRefList>
        </Allele>
    """)

    result = ClinvarVariationAnnotator._extract_xrefs(soup)

    assert result['uniprot_id'] == 'Uniprot-1'
    assert result['omim_id'] == 'Omim-1'
    assert result['dbsnp_id'] == 'rs123'


def test_extract_molecular_consequences():
    soup = make_soup("""
        <Allele>
            <MolecularConsequenceList>
                <MolecularConsequence HGVS="change-1"
                                      Function="type-1" />
            </MolecularConsequenceList>
        </Allele>
    """)

    results = ClinvarVariationAnnotator._extract_molecular_consequences(soup)

    result = results[0]
    assert result['hgvs'] == 'change-1'
    assert result['function'] == 'type-1'


def test_extract_allele_frequencies():
    soup = make_soup("""
        <Allele>
            <AlleleFrequency Value="0.01" Type="Type-1" MinorAllele="A" />
            <AlleleFrequency Value="0.02" Type="Type-2" MinorAllele="A" />
        </Allele>
    """)

    result = ClinvarVariationAnnotator._extract_allele_frequencies(soup)

    assert result == {'A': {'Type-1': 0.01, 'Type-2': 0.02}}


def test_generate_clinical_summary():
    assertions = [
        {'clinical_significance': 'ClinSig-1'},
        {'clinical_significance': 'ClinSig-1'},
        {'clinical_significance': 'ClinSig-2'},
    ]

    result = ClinvarVariationAnnotator._generate_clinical_summary(assertions)
    assert result == {'ClinSig-1': 2, 'ClinSig-2': 1}

