import pandas as pd
import pytest

from anotala.pipeline import annotate_clinvar_accessions


@pytest.mark.skip('Not implemented completely')
def test_annotate_clinvar_accessions():
    clinvar_entries = pd.Series([
        {'accession': 'ACC-1'},
        {'accession': 'ACC-2'},
    ])

    annotate_clinvar_accessions(clinvar_entries)

    # Stub ClinvarRCVAnnotator
    # check call to .annotate(['ACC-1', 'ACC-2'])
    # stub return-value
    # check that value is added to the entries
    #  assert 0

