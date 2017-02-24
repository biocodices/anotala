import pytest

from anotamela.helpers import (
    read_VEP_header,
    read_VEP_table,
)


@pytest.fixture
def tsv():
    return pytest.helpers.file('vep.tsv')


def test_read_VEP_header(tsv):
    fields = read_VEP_header(tsv)
    expected_fields = 'Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	IMPACT	DISTANCE	STRAND	FLAGS	VARIANT_CLASS	SYMBOL	SYMBOL_SOURCE	HGNC_ID	BIOTYPE	CANONICAL	TSL	APPRIS	CCDS	ENSP	SWISSPROT	TREMBL	UNIPARC	REFSEQ_MATCH	GENE_PHENO	SIFT	PolyPhen	EXON	INTRON	DOMAINS	HGVSc	HGVSp	HGVS_OFFSET	GMAF	AFR_MAF	AMR_MAF	EAS_MAF	EUR_MAF	SAS_MAF	AA_MAF	EA_MAF	ExAC_MAF	ExAC_Adj_MAF	ExAC_AFR_MAF	ExAC_AMR_MAF	ExAC_EAS_MAF	ExAC_FIN_MAF	ExAC_NFE_MAF	ExAC_OTH_MAF	ExAC_SAS_MAF	CLIN_SIG	SOMATIC	PHENO	PUBMED	MOTIF_NAME	MOTIF_POS	HIGH_INF_POS	MOTIF_SCORE_CHANGE'.split()
    assert fields == expected_fields


def test_read_VEP_table(tsv):
    df = read_VEP_table(tsv, keep_only_canonical=True)
    assert len(df) == 2
    assert list(df['CANONICAL']) == ['YES', 'YES']

    df = read_VEP_table(tsv, keep_only_canonical=False)
    assert len(df) == 16

