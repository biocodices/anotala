import pytest

from anotamela.helpers import rsids_from_vcf


@pytest.fixture
def vcf_path():
    return pytest.helpers.file('1_sample.vcf')


def test_rsids_from_vcf(vcf_path):
    result = rsids_from_vcf(vcf_path)
    assert len(result) == 1
    assert result == ['rs268']
