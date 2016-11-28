import pytest

from anotamela.annotators import AnnotatorWithCache


def test_annotate():
    ids_to_annotate = '1 2 3'.split()

    generic_annotator = AnnotatorWithCache(cache='redis')
    # WIP!

