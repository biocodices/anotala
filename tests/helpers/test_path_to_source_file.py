import pytest

from anotamela.helpers import path_to_source_file


def test_path_to_source_file():
    result = path_to_source_file('clinvar_20180729.grch37.vcf.gz')
    assert result.endswith('anotamela/anotamela/sources/clinvar_20180729.grch37.vcf.gz')

    with pytest.raises(Exception):
        path_to_source_file('non-existent.file')
