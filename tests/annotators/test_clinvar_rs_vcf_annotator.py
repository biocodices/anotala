import pandas as pd
import pytest

from anotamela.annotators import ClinvarRsVCFAnnotator


# With the real ClinVar VCF, this test would take too long
def path_to_test_clinvar_vcf():
    return pytest.helpers.file('clinvar.vcf')


def test_init():
    annotator = ClinvarRsVCFAnnotator(path_to_test_clinvar_vcf())
    assert len(annotator.data) == 3


def test_filter_data_by():
    df = pd.DataFrame([
        {'rs_id': 'rs1', 'chrom': '1'},
        {'rs_id': 'rs1', 'chrom': '1'},
        {'rs_id': 'rs2', 'chrom': '2'},
        {'rs_id': 'rs3', 'chrom': '3'},
    ])
    result = ClinvarRsVCFAnnotator._filter_data_by(df, 'rs_id', ['rs1', 'rs3'])
    assert len(result) == 3
    assert result['rs_id'].values.tolist() == ['rs1', 'rs1', 'rs3']


def test_annotate():
    ids_to_annotate = ['rs1', 'rs2']
    annotator = ClinvarRsVCFAnnotator(path_to_test_clinvar_vcf())
    result = annotator.annotate(ids_to_annotate)
    assert len(result['rs1']) == 2
    assert [record['variation_id'] for record in result['rs1']] == ['Var-1', 'Var-2']
    assert len(result['rs2']) == 1
