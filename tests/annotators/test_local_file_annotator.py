from unittest.mock import Mock
import pytest

from anotamela.annotators.base_classes import LocalFileAnnotator


@pytest.fixture
def annotator():
    return LocalFileAnnotator(path_to_annotations_file='/foo/bar')


def test_init(annotator):
    assert annotator.path == '/foo/bar'


def test_data(annotator):
    mock_load_data = Mock()
    annotator._load_data = mock_load_data

    # This should call _load_data the first time:
    annotator.data
    mock_load_data.assert_called_once_with()

    # Data should be stored in an internal instance variable:
    assert hasattr(annotator, '_data')

    # This should not call _load_data the second time:
    annotator.data
    mock_load_data.assert_called_once_with()
