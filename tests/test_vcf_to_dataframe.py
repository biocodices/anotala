from os.path import dirname, realpath, join
import re
import gzip

import pytest
from pytest import raises

from anotamela import vcf_to_dataframe


TEST_DIR = dirname(realpath(__file__))


def _test_file(filename):
    return join(TEST_DIR, filename)

def _test_vcf_to_dataframe(filename, gzipped):
    vcf_path = _test_file(filename)

    df = vcf_to_dataframe(vcf_path)
    assert df.shape == (1742, 9)
    assert len(df['chrom'].cat.categories) == 5

    # Read the file "manually" to compare the data
    if gzipped:
        with gzip.open(vcf_path) as f:
            rows = [line.decode('utf-8') for line in f.readlines()
                    if line and not line.startswith(b'#')]
    else:
        with open(vcf_path) as f:
            rows = [line for line in f.readlines()
                    if line and not line.startswith('#')]

    records = [re.split(r'\s+', row) for row in rows]
    seen_ids = set(record[2] for record in records)

    assert seen_ids == set(df['id'].unique())

    df = vcf_to_dataframe(vcf_path, keep_samples='HG00096')
    assert 'HG00096' in df.columns
    assert df.shape == (1742, 10)

    with raises(ValueError):
        vcf_to_dataframe(vcf_path, keep_samples='non existent')

def test_vcf_to_dataframe_gzipped():
    _test_vcf_to_dataframe('files/sample.vcf.gz', gzipped=True)

def test_vcf_to_dataframe():
    _test_vcf_to_dataframe('files/sample.vcf', gzipped=False)

