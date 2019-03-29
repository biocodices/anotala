from anotala.annotators.base_classes import Annotator


def test_set_of_string_ids():
    ids = [1, '2']
    result = Annotator._set_of_string_ids(ids)
    assert result == {'1', '2'}
