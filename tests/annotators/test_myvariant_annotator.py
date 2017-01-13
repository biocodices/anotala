from anotamela.annotators.base_classes import MyVariantAnnotator


def test_myvariant_annotator_parse_annotationa(monkeypatch):
    # Test the results are returned in the same order
    hits = [
        {'_id': 'chr1:123A>G', 'expected_ix': 1},
        {'_id': 'chr1:123A>C', 'expected_ix': 0},
        {'_id': 'chr1:123A>T', 'expected_ix': 2},
    ]

    annotations = MyVariantAnnotator._parse_annotation(hits)
    assert [r['expected_ix'] for r in annotations] == [0, 1, 2]

    # If each of the parsed hits is already a list, the parsing should merge
    # them into a single flat list:

    # We mock the classmethod to parse hits to get the list of lists
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

