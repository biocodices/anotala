import pytest

from anotamela.annotators import DbSNPAnnotator


def test_dbsnp_annotator():
    ids_to_annotate = 'rs268 rs123'.split()

    annotator = DbSNPAnnotator(cache='mock_cache')
    info_dict = annotator.annotate(ids=ids_to_annotate)

    assert all(id_ in info_dict for id_ in ids_to_annotate)
    assert all(response for response in info_dict.values())  # No empty responses

