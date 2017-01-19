import pytest

from anotamela.pipeline import read_variants_from_vcf


def test_read_variants_from_vcf():
    vcf = pytest.helpers.file('test.vcf')
    result = read_variants_from_vcf(vcf)

    assert len(result['rs_variants']) == 1
    assert len(result['other_variants']) == 1

