from anotamela.annotators.base_classes import MyVariantAnnotator


def test_myvariant_annotator_parse_annotation():
    # Test the results are returned in the same order
    hits = [
        {'_id': 'chr1:123A>G', 'expected_ix': 1},
        {'_id': 'chr1:123A>C', 'expected_ix': 0},
        {'_id': 'chr1:123A>T', 'expected_ix': 2},
    ]

    annotator = MyVariantAnnotator('mock_cache')
    results = annotator._parse_annotation(hits)
    assert [r['expected_ix'] for r in results] == [0, 1, 2]

