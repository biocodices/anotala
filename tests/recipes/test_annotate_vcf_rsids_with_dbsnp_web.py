import os
import json
import pytest

from anotamela.recipes.annotate_vcf_rsids_with_dbsnp_web import (
    annotate_rsids_with_dbsnp_web,
    annotate_vcf_rsids_with_dbsnp_web
)


@pytest.fixture
def vcf_path():
    return pytest.helpers.file('1_sample.vcf')


def test_annotate_vcf_rsids_with_dbsnp_web(vcf_path, dict_cache, proxies):
    result = annotate_vcf_rsids_with_dbsnp_web(
        vcf_path,
        cache=dict_cache,
        proxies=proxies
    )
    assert result[0]['rs_id'] == 'rs268'

    out_fn = '/tmp/__anotamela_test_annotate_vcf_rsids_with_dbsnp_web.json'
    try:
        result = annotate_vcf_rsids_with_dbsnp_web(
            vcf_path,
            output_json_path=out_fn,
            cache=dict_cache
        )
        with open(out_fn) as f:
            annotations = json.load(f)

        assert annotations[0]['rs_id'] == 'rs268'
    finally:
        os.remove(out_fn)


def test_annotate_rsids_with_dbsnp_web(dict_cache):
    result = annotate_rsids_with_dbsnp_web(['rs268'], cache=dict_cache)
    assert result[0]['rs_id'] == 'rs268'
