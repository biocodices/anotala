from os.path import join
import pickle

import pytest

from anotamela.annotators import *
from helpers import get_test_file


ANNOTATOR_CLASSES = [
    ClinvarRsAnnotator,
    DbsnpMyvariantAnnotator,
    DbsnpEntrezAnnotator,
    DbsnpWebAnnotator,
    HgvsAnnotator,
    GeneEntrezAnnotator,
    MafAnnotator,
    MygeneAnnotator,
    OmimGeneAnnotator,
    OmimVariantAnnotator,
    PubmedAnnotator,
    SnpeffAnnotator,
    UniprotAnnotator,
]


@pytest.mark.parametrize('annotator_class', ANNOTATOR_CLASSES)
def test_parse_annotation(annotator_class):
    if not hasattr(annotator_class, '_parse_annotation'):
        return

    name = annotator_class.SOURCE_NAME
    raw_annotation = get_annotation(name, parsed=False)
    parsed_annnotation_expected = get_annotation(name, parsed=True)
    parsed_annnotation_seen = \
        annotator_class._parse_annotation(raw_annotation)

    assert parsed_annnotation_seen == parsed_annnotation_expected


def get_annotation(source_name, parsed):
    subdir = 'parsed_annotations' if parsed else 'raw_annotations'
    filename = join(subdir, source_name + '.pickle')
    with open(get_test_file(filename), 'rb') as f:
        return pickle.load(f)

