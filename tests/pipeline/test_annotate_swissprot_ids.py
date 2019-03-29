from anotala.annotators import UniprotAnnotator
from anotala.pipeline import annotate_swissprot_ids


def test_annotate_swissprot_ids(monkeypatch):
    def mock_annotate(self, ids, **kwargs):
        # Returns 2 fake variants per ID
        return {id_: [{'foo': 'bar'}, {'baz': 'qux'}]
                for id_ in ids}

    monkeypatch.setattr(UniprotAnnotator, 'annotate', mock_annotate)

    ids = ['SP1', 'SP2']
    result = annotate_swissprot_ids(ids, cache='mock_cache', use_web=True,
                                    use_cache=True, proxies={'foo': 'bar'})

    assert len(result) == 4  # 2 fake variants * 2 ids
    assert {'foo': 'bar'} in result
    assert {'baz': 'qux'} in result

