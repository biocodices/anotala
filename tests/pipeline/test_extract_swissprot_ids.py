import numpy as np
import pandas as pd

from anotala.pipeline import extract_swissprot_ids


def test_extract_swissprot_ids():
    assert extract_swissprot_ids(pd.Series([])) == []

    mygene_annotations = pd.Series([
        None,
        np.nan,
        {},  # No 'swissprot' in this gene entry
        {'swissprot': 'SP1'},  # Single ID
        {'swissprot': ['SP2']},  # List of IDs
        {'swissprot': ['SP1', 'SP3']},
    ])

    swissprot_ids = extract_swissprot_ids(mygene_annotations)
    assert sorted(swissprot_ids) == ['SP1', 'SP2', 'SP3']

