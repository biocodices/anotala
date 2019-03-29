from unittest.mock import Mock
import pytest

from anotala.annotators.base_classes import LocalFileAnnotator


@pytest.fixture
def annotator():
    return LocalFileAnnotator(path_to_annotations_file='/foo/bar')


def test_init(annotator):
    assert annotator.path == '/foo/bar'

    class ChildAnnotator(LocalFileAnnotator):
        PATH_TO_ANNOTATIONS_FILE = "/bar/baz"

    custom_file_local_annotator = ChildAnnotator()
    assert custom_file_local_annotator.path == "/bar/baz"


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
