from bs4 import BeautifulSoup

import anotamela
from anotamela.annotators import ClinvarVariationAnnotator


# See https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&id=1550&rettype=variation
# For an example of the kind of XML structure we are dealing with.

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
    soup = make_soup('<VariationReport VariationID="Var-1"></VariationReport>')
    result = ClinvarVariationAnnotator._extract_variation_id(soup)
    assert result == 'Var-1'


def test_extract_variation_name():
    soup = make_soup('<VariationReport VariationName="Var-1-Name"></VariationReport>')
    result = ClinvarVariationAnnotator._extract_variation_name(soup)
    assert result == 'Var-1-Name'


def test_extract_variation_type():
    soup = make_soup('<VariationReport VariationType="Var-1-Type"></VariationReport>')
    result = ClinvarVariationAnnotator._extract_variation_type(soup)
    assert result == 'Var-1-Type'

    soup = make_soup('<VariationReport></VariationReport>')
    # Just test it does not break when VariantType is missing
    assert ClinvarVariationAnnotator._extract_variation_type(soup) is None


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


def test_extract_single_gene_name():
    soup = make_soup("""
        <VariationReport>
            <GeneList GeneCount="1">
                <Gene Symbol="Gene-Symbol-1"></Gene>
            <GeneList/>
        </VariationReport>
    """)

    result = ClinvarVariationAnnotator._extract_single_gene_name(soup)
    assert result == 'Gene-Symbol-1'

    soup = make_soup("""
        <VariationReport>
            <GeneList GeneCount="2">
                <Gene Symbol="Gene-Symbol-1"></Gene>
                <Gene Symbol="Gene-Symbol-2"></Gene>
            <GeneList/>
        </VariationReport>
    """)

    result = ClinvarVariationAnnotator._extract_single_gene_name(soup)
    assert result is None

    soup = make_soup("""<VariationReport></VariationReport>""")
    result = ClinvarVariationAnnotator._extract_single_gene_name(soup)
    assert result is None


def test_parse_phenotype_list(monkeypatch):
    phenotype_list = make_soup("""
        <PhenotypeList>
            <Phenotype Name="Pheno-1"></Phenotype>
            <Phenotype Name="Pheno-2">
                <XRefList>
                    <XRef ID="MIM-2" DB="OMIM" />
                </XRefList>
            </Phenotype>
        </PhenotypeList>
    """)

    result = ClinvarVariationAnnotator._parse_phenotype_list(phenotype_list)
    assert len(result) == 2
    assert result[0] == {'name': 'Pheno-1'}
    assert result[1] == {'name': 'Pheno-2', 'omim_id': 'MIM-2',
                         'incidental': False}


def test_extract_obervation():
    soup = make_soup("""
        <VariationReport VariationID="1">
            <ObservationList>
                <Observation VariationID="1" ObservationType="primary">
                    <ClinicalSignificance DateLastEvaluated="2001-01-01">
                        <Description>Pathogenic</Description>
                    </ClinicalSignificance>
                    <ReviewStatus>criteria provided</ReviewStatus>
                    <PhenotypeList>
                        <Phenotype Name="Pheno-1"></Phenotype>
                    </PhenotypeList>
                </Observation>
                <Observation VariationID="2">
                    <ClinicalSignificance>
                        <Description>Benign</Description>
                    </ClinicalSignificance>
                </Observation>
            </ObservationList>
        </VariationReport>
    """)
    result = ClinvarVariationAnnotator._extract_observation(soup)
    assert result == {
        'variation_id': '1',
        'clinical_significance': 'Pathogenic',
        'date_last_evaluated': '2001-01-01',
        'type': 'primary',
        'review_status': 'criteria provided',
        'phenotypes': [{'name': 'Pheno-1'}],
    }


def test_extract_clinical_assertions():
    soup = make_soup("""
        <VariationReport>
            <ClinicalAssertionList>
                <GermlineList>
                    <Germline SubmitterName="Submitter-1"
                              DateLastSubmitted="2000-01-01">

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
                <SomaticList>
                    <Somatic SubmitterName="Submitter-2"
                             DateLastSubmitted="2000-01-02">

                        <PhenotypeList>
                            <Phenotype Name="Pheno-2">
                                <XRefList>
                                    <XRef Type="MIM" ID="MIM-2" DB="OMIM"/>
                                </XRefList>
                            </Phenotype>
                        </PhenotypeList>

                        <ClinicalSignificance>
                            <Description>ClinSig-2</Description>
                            <Method>Method-2</Method>
                        </ClinicalSignificance>

                    </Somatic>
                </SomaticList>
            </ClinicalAssertionList>
        </VariationReport>
    """)

    results = ClinvarVariationAnnotator._extract_clinical_assertions(soup)
    assert len(results) == 2

    germline = results[0]
    assert germline == {
        'clinical_significance': 'ClinSig-1',
        'date_last_submitted': '2000-01-01',
        'method': 'Method-1',
        'phenotypes': [{'name': 'Pheno-1', 'omim_id': 'MIM-1',
                        'incidental': False}],
        'submitter_name': 'Submitter-1',
        'type': 'germline',
    }

    somatic = results[1]
    assert somatic == {
        'clinical_significance': 'ClinSig-2',
        'date_last_submitted': '2000-01-02',
        'method': 'Method-2',
        'phenotypes': [{'name': 'Pheno-2', 'omim_id': 'MIM-2',
                        'incidental': False}],
        'submitter_name': 'Submitter-2',
        'type': 'somatic',
    }


    soup = make_soup('<VariationReport><GermlineList><Germline></Germline>'
                     '</GermlineList></VariationReport>')
    # Make sure it does not break on missing data
    ClinvarVariationAnnotator._extract_clinical_assertions(soup)


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

    soup = make_soup("""
        <Allele>
            <SequenceLocation Assembly="GRCh37"
                              innerStart="1000"
                              innerStop="1000" />
        </Allele>
    """)
    result = ClinvarVariationAnnotator._extract_sequence_info_from_allele(soup)
    assert result['start_g37'] == 1000
    assert result['stop_g37'] == 1000

def test_extract_allele_hgvs():
    soup = make_soup("""
        <Allele>
            <HGVSlist>
                <HGVS Assembly="GRCh37"
                      Change="genomic-change-37"
                      AccessionVersion="Accession-1">Name-37</HGVS>
                <HGVS Assembly="GRCh38"
                      Change="genomic-change-38"
                      AccessionVersion="Accession-1">Name-38</HGVS>
                <HGVS Type="HGVS, coding, RefSeq">accession:coding-change-1</HGVS>
                <HGVS Type="HGVS, coding, RefSeq">accession:coding-change-2</HGVS>
                <HGVS Type="HGVS, protein, RefSeq">accession:protein-change-1</HGVS>
            </HGVSlist>
        </Allele>
    """)

    result = ClinvarVariationAnnotator._extract_allele_hgvs(soup)

    assert result['genomic_change_g37'] == 'genomic-change-37'
    assert result['genomic_change_g37_accession'] == 'Accession-1'
    assert result['genomic_change_g37_name'] == 'Name-37'

    assert result['genomic_change_g38'] == 'genomic-change-38'
    assert result['genomic_change_g38_accession'] == 'Accession-1'
    assert result['genomic_change_g38_name'] == 'Name-38'

    assert result['coding_changes'] == ['accession:coding-change-1',
                                        'accession:coding-change-2']

    assert result['protein_changes'] == ['accession:protein-change-1']

    soup = make_soup("<Allele></Allele>")
    # Test it does not break when data is missing
    ClinvarVariationAnnotator._extract_allele_hgvs(soup)


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

    soup = make_soup("""
        <Allele>
            <XRefList>
            </XRefList>
        </Allele>
    """)
    # Check it does not break when info is missing
    ClinvarVariationAnnotator._extract_xrefs(soup)


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


def test_associated_phenotypes():
    clinical_assertions = [
        {'phenotypes': [{'name': 'Pheno-1'}, {'name': 'Pheno-2'}]},
        {'phenotypes': [{'name': 'Pheno-2'}, {'name': 'Pheno-3'}]},
        {} # Has no phenotypes key, shouldn't break
    ]

    result = ClinvarVariationAnnotator._associated_phenotypes(clinical_assertions)
    assert result == ['Pheno-1', 'Pheno-2', 'Pheno-3']

