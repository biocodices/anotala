import pandas as pd

from anotamela.annotators import ClinvarRsVCFAnnotator


def test_init(path_to_test_clinvar_vcf):
    annotator = ClinvarRsVCFAnnotator(path_to_test_clinvar_vcf)
    assert len(annotator.data) == 4


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


def test_annotate(path_to_test_clinvar_vcf):
    ids_to_annotate = ['rs1', 'rs2']
    annotator = ClinvarRsVCFAnnotator(path_to_test_clinvar_vcf)
    result = annotator.annotate(ids_to_annotate)
    assert len(result['rs1']) == 1
    assert [record['variation_id'] for record in result['rs1']] == ['Var-1']
    assert len(result['rs2']) == 1
