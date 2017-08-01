import numpy as np

from anotamela.pipeline import extract_ensembl_consequence


def test_extract_ensembl_consequence():
    result = extract_ensembl_consequence({'most_severe_consequence': 'foo'})
    assert result == 'foo'

    result = extract_ensembl_consequence({'most_severe_consequence': np.nan})
    assert result == None

    result = extract_ensembl_consequence({'foo': 'bar'})
    assert result == None

    result = extract_ensembl_consequence(None)
    assert result == None

