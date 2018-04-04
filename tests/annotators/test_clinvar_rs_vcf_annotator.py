import pytest

from anotamela.annotators import ClinvarRsVCFAnnotator


def path_to_vcf():
    return pytest.helpers.file('clinvar.vcf')


def test_init():
    annotator = ClinvarRsVCFAnnotator(path_to_vcf())
    assert len(annotator.data) == 3

    # This test takes too long, the VCF is big!
    # With package-provided ClinVar VCF:
    #  annotator = ClinvarRsVCFAnnotator()
    #  assert len(annotator.data) == 341_303


def test_annotate():
    ids_to_annotate = ['rs1', 'rs2']
    annotator = ClinvarRsVCFAnnotator(path_to_vcf())
    result = annotator.annotate(ids_to_annotate)
    assert len(result['rs1']) == 2
    assert [record['variation_id'] for record in result['rs1']] == ['Var-1', 'Var-2']
    assert len(result['rs2']) == 1
