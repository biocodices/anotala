from anotamela.annotators.base_classes import AnnotatorWithCache


def test_annotator_with_cache():
    annotator = AnnotatorWithCache('mock_cache')
    inexistent_id = 'inexistent_id'
    assert annotator.annotate_one(inexistent_id, use_web=False) is None
    assert annotator.annotate(inexistent_id, use_web=False) == {}

