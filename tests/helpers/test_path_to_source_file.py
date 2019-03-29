import pytest

from anotala.helpers import path_to_source_file


def test_path_to_source_file():
    fn = 'integrated_call_samples_v3.20130502.ALL.panel'
    result = path_to_source_file(fn)
    assert result.endswith(f'anotala/anotala/sources/{fn}')

    with pytest.raises(Exception):
        path_to_source_file('non-existent.file')
