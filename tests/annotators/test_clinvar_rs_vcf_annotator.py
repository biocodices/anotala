import pandas as pd
import pytest

from anotamela.annotators import ClinvarRsVCFAnnotator


def path_to_vcf():
    return pytest.helpers.file('clinvar.vcf')


def test_init():
    annotator = ClinvarRsVCFAnnotator(path_to_vcf())
    assert len(annotator.data) == 2


def test_read_vcf():
    result = ClinvarRsVCFAnnotator._read_vcf(path_to_vcf())
    assert len(result) == 2


def test_parse_dataframe():
    df = pd.DataFrame([
        {'info': {'RS': 'rs1',
                  'CLNHGVS': 'name-1',
                  'FOO': 'foo'}},
        {'info': {'RS': 'rs2',
                  'CLNHGVS': 'name-2',
                  'BAR': 'bar'}},
        {'info': {'RS': 'rs3',
                  'CLNHGVS': 'name-3'}},
    ])
    ClinvarRsVCFAnnotator._parse_dataframe(df)
    assert df['rs_id'].values.tolist() == ['rs1', 'rs2', 'rs3']
    assert df['clinvar_name'].values.tolist() == ['name-1', 'name-2', 'name-3']
    assert df['FOO'].values.tolist() == ['foo', None, None]
    assert df['BAR'].values.tolist() == [None, 'bar', None]


def test_annotate_many_ids():
    ids_to_annotate = ['rs1', 'rs2']
    annotator = ClinvarRsVCFAnnotator(path_to_vcf())
    result = annotator._annotate_many_ids(ids_to_annotate)
    assert result['rs1'] == {'clinvar_variation_id': 'Clinvar-1'}
    assert result['rs2'] == {'clinvar_variation_id': 'Clinvar-2'}
