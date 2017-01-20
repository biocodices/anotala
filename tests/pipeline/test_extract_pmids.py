from anotamela.pipeline import extract_pmids


def test_extract_pmids():
    omim_entries = [
        {'pubmed_entries': [{'pmid': 1}, {'pmid': 2}]},
        {'pubmed_entries': [{'pmid': 1}, {'pmid': 3}]},
        {'pubmed_entries': []},
        {},  # Doesn't have a 'pubmed_entries' key
    ]

    pmids = extract_pmids(omim_entries)
    assert pmids == [1, 2, 3]

