import pytest

from anotamela.annotators import (
        DbsnpWebAnnotator,
        DbsnpEntrezAnnotator,
    )


ids_to_annotate = 'rs268 rs123'.split()

@pytest.mark.parametrize('annotator_class', [
        DbsnpWebAnnotator,
        DbsnpEntrezAnnotator
    ])

def test_generic_annotator(annotator_class):
    annotator = annotator_class(cache='mock_cache')
    info_dict = annotator.annotate(ids=ids_to_annotate)

    assert all(id_ in info_dict for id_ in ids_to_annotate)
    assert all(response for response in info_dict.values())  # No empty responses

