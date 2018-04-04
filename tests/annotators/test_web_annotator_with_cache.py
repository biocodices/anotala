from anotamela.annotators.base_classes import WebAnnotatorWithCache


def test_annotator_with_cache():
    annotator = WebAnnotatorWithCache('mock_cache')
    inexistent_id = 'inexistent_id'
    assert annotator.annotate_one(inexistent_id, use_web=False) is None
    assert annotator.annotate(inexistent_id, use_web=False) == {}

