import pandas as pd

from anotamela.annotators import ClinvarPosVCFAnnotator


def test_filter_data_by():
    df = pd.DataFrame([
        {'rs_id': 'rs1', 'chrom_pos': '1:100',
         'variant_type': 'single_nucleotide_variant'},
        {'rs_id': 'rs1a', 'chrom_pos': '1:100',
         'variant_type': 'single_nucleotide_variant'},
        {'rs_id': 'rs1b', 'chrom_pos': '1:110',
         'variant_type': 'Deletion'},
        {'rs_id': 'rs2', 'chrom_pos': '2:200',
         'variant_type': 'single_nucleotide_variant'},
        {'rs_id': 'rs3', 'chrom_pos': '3:300',
         'variant_type': 'single_nucleotide_variant'},
    ])

    positions = ['1:100', '1:110', '3:300']
    result = ClinvarPosVCFAnnotator._filter_data_by(df, 'chrom_pos', positions)
    assert len(result) == 4
    assert result['rs_id'].values.tolist() == ['rs1', 'rs1a', 'rs1b', 'rs3']


def test_annotate(path_to_test_clinvar_vcf):
    positions_to_annotate = ['1:100', '1:110', '2:200']
    annotator = ClinvarPosVCFAnnotator(path_to_test_clinvar_vcf)
    result = annotator.annotate(positions_to_annotate)
    assert len(result['1:100']) == 1
    assert len(result['1:110']) == 1
    assert len(result['2:200']) == 1
    assert [record['variation_id'] for record in result['1:100']] == ['Var-1']
