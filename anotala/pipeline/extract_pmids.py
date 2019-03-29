def extract_pmids(entries):
    """
    Given a list of entries with 'pubmed_entries' inside, extact the unique
    PubMed IDs mentioned.
    """
    pmids = set()

    for entry in entries:
        for pubmed in entry.get('pubmed_entries', []):
            pmids.add(pubmed['pmid'])

    return list(pmids)

