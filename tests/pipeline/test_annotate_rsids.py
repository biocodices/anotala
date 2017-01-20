from anotamela.annotators import RS_ANNOTATOR_CLASSES
from anotamela.pipeline import annotate_rsids


def test_annotate_rsids(monkeypatch):
    def mock_annotate(_, ids, *args, **kwargs):
        fake_annotations = ['{}-ann'.format(id_) for id_ in ids]
        return dict(zip(ids, fake_annotations))

    annotator_classes = RS_ANNOTATOR_CLASSES.values()

    for annotator_class in annotator_classes:
        monkeypatch.setattr(annotator_class, 'annotate', mock_annotate)

    ids = ['rs1', 'rs2', 'rs3']
    result = annotate_rsids(ids, cache='mock_cache',
                            use_web=True, use_cache=True,
                            annotator_names=RS_ANNOTATOR_CLASSES.keys(),
                            proxies={'foo': 'bar'})

    # Check that the result dataframe is populated with all the annotators data
    for annotator_class in annotator_classes:
        assert annotator_class.SOURCE_NAME in result

