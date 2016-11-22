from os.path import dirname, realpath, join
import re

import pytest

from anotamela import variants_from_vcf


TEST_DIR = dirname(realpath(__file__))


def _test_file(filename):
    return join(TEST_DIR, filename)


def test_variants_from_vcf():
    vcf_path = _test_file('files/sample.vcf')

    with open(vcf_path) as f:
        variants = [line for line in f.read().split('\n')
                    if line and not line.startswith('#')]

    variants = [re.split(r'\s+', variant) for variant in variants]
    seen_ids = set(variant[2] for variant in variants)
    assert seen_ids == set(variants_from_vcf(vcf_path))

