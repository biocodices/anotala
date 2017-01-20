def extract_pmids(omim_entries):
    """
    Given a list of omim variant entries, extact the unique PubMed IDs cited.
    """
    pmids = set()

    for entry in omim_entries:
        for pubmed in entry.get('pubmed_entries', []):
            pmids.add(pubmed['pmid'])

    return list(pmids)

