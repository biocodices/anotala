from copy import deepcopy


def assign_pubmed_entries_to_omim_entries(omim_entries, pubmed_entries_by_pmid):
    """
    Given a list of *omim_entries*, which have 'pubmed_entries' in them,
    update the info of each PubMed entry with the annotations by PMID in
    passed in *pubmed_entries_by_pmid*.

    Return a new list of omim_entries with the extra data added.
    """
    updated_omim_entries = deepcopy(omim_entries)

    for entry in updated_omim_entries:
        for pubmed_entry in entry.get('pubmed_entries', []):
            annotation = pubmed_entries_by_pmid.get(pubmed_entry['pmid'])
            if annotation:
                pubmed_entry.update(annotation)

    return updated_omim_entries

