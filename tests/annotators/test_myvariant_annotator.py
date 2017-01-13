from anotamela.annotators.base_classes import MyVariantAnnotator


def test_myvariant_annotator_parse_annotationa(monkeypatch):
    # Test the results are returned in the expected order (by _id)
    hits = [
        {'_id': 'B', 'expected_ix': 1},
        {'_id': 'A', 'expected_ix': 0},
        {'_id': 'C', 'expected_ix': 2},
    ]

    annotations = MyVariantAnnotator._parse_annotation(hits)
    assert [r['expected_ix'] for r in annotations] == [0, 1, 2]

    # If all the parsed hits are lists, the parsing should merge them into a
    # single flat list

    # We mock the classmethod _parse_hit to get the list inside annotations
    monkeypatch.setattr(MyVariantAnnotator, '_parse_hit',
                        lambda hit: hit['test-list'])

    hits = [{'_id': 'A',
             'test-list': [{'foo': 'bar'},
                           {'bar': 'baz'}]},
            {'_id': 'B',
             'test-list': [{'qux': 'boo'}]}]

    # Here the mocked class method will be callled for each of the hits:
    annotations = MyVariantAnnotator._parse_annotation(hits)

    # The three annotations inside 'test-list' keys should be in a flat list:
    assert len(annotations) == 3

