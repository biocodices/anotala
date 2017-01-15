import pytest

from anotamela import SnpeffAnnotator


def test_snpeff_annotator_parse_hit():
    hit_single = {'snpeff': {'ann': {'effect': 'eff-1&eff-2'}}}
    annotations = SnpeffAnnotator._parse_hit(hit_single)

    assert isinstance(annotations, list)
    assert len(annotations) == 1
    assert 'effect' not in annotations[0]

    hit_multiple = {'snpeff': {'ann': [{}, {}, {}]}}
    annotations = SnpeffAnnotator._parse_hit(hit_multiple)
    assert len(annotations) == 3


@pytest.mark.parametrize('annotation,expected_effects', [
    ({'effect': 'eff-1'}, ['eff-1']),
    ({'effect': 'eff-1&eff-2'}, ['eff-1', 'eff-2']),
    ({}, []),
])
def test_split_effects(annotation, expected_effects):
    assert SnpeffAnnotator._split_effects(annotation) == expected_effects


@pytest.mark.parametrize('annotation,expected_allele', [
    ({'hgvs_c': 'NM_1.1:c.123A>T'}, 'T'),
    ({'hgvs_c': 'NM_015001.2:c.*11199_*11199delTTT'}, 'delTTT'),
    ({}, None),
])
def test_infer_allele(annotation, expected_allele):
    assert SnpeffAnnotator._infer_allele(annotation) == expected_allele

