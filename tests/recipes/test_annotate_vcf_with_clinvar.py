import pytest
import json
import os

from anotamela.recipes.annotate_vcf_rsids_with_clinvar import (
    rsids_from_vcf,
    annotate_rsids_with_clinvar,
    annotate_vcf_rsids_with_clinvar
)


@pytest.fixture
def vcf_path():
    return pytest.helpers.file('1_sample.vcf')


def test_annotate_vcf_rsids_with_clinvar(vcf_path, dict_cache):
    result = annotate_vcf_rsids_with_clinvar(
        vcf_path,
        cache=dict_cache
    )
    assert result[0]['dbsnp_id'] == 'rs268'

    out_fn = '/tmp/__anotamela_test_annotate_vcf_rsids_with_clinvar.json'
    try:
        result = annotate_vcf_rsids_with_clinvar(
            vcf_path,
            output_json_path=out_fn,
            cache=dict_cache
        )
        with open(out_fn) as f:
            annotations = json.load(f)

        assert annotations[0]['dbsnp_id'] == 'rs268'
    finally:
        os.remove(out_fn)


def test_annotate_rsids_with_clinvar(dict_cache):
    result = annotate_rsids_with_clinvar(['rs268'], cache=dict_cache)
    assert result[0]['dbsnp_id'] == 'rs268'


def test_rsids_from_vcf(vcf_path):
    result = rsids_from_vcf(vcf_path)
    assert len(result) == 1
    assert result == ['rs268']
