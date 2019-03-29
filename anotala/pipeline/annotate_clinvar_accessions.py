from copy import deepcopy
from operator import itemgetter

import pandas as pd

from anotala.annotators import ClinvarRCVAnnotator


def annotate_clinvar_accessions(clinvar_entries):
    """
    Expects a series 'clinvar_entries'. Each value would be a list of
    dictionaries.

    First, it annotates all the 'accession' RCV ids found across all entries.

    Then, it adds the accession-related data to the original entries.

    Returns a copy of the series, with the extra RCV accession info added to
    each entry.
    """

    # WIP CODE!

    #  new_series = pd.Series([])

    #  accessions = clinvar_entries.map(itemgetter('accession'))
    #  accessions_annotations = ClinvarRCVAnnotator.annotate(accessions)

    #  for ix, entries in clinvar_entries.items():
        #  new_entries = deepcopy(entries)

        #  # Add the new RCV annotations to the entries
        #  for new_entry in new_entries:
            #  accession = new_entry.get('accession')
            #  annotation = accessions_annotations.get(accession, {})
            #  new_entry.update(annotation)

        #  #  accession_data = accessions_annotations[]
        #  new_entries.update()
        #  # new_series[ix] =
