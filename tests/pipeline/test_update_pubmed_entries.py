from anotala.pipeline import update_pubmed_entries


def test_update_pubmed_entries():
    omim_entries = [
        {'id': 1},  # Has no 'pubmed_entries'
        {'id': 2, 'pubmed_entries': []},
        {'id': 3, 'pubmed_entries': [{'pmid': 1}, {'pmid': 2}]},
        {'id': 4, 'pubmed_entries': [{'pmid': 1}, {'pmid': 3}]},
    ]

    pubmed_entries = {
        1: {'foo': 'bar'},
        2: {'foo': 'baz'},
        # No annotation for pmid 3
    }

    new_omim_entries = update_pubmed_entries(omim_entries, pubmed_entries)

    for omim_entry in new_omim_entries:
        id_ = omim_entry['id']

        if id_ == 1:
            assert 'pubmed_entries' not in omim_entry
        if id_ == 2:
            assert omim_entry['pubmed_entries'] == []
        if id_ == 3:
            assert omim_entry['pubmed_entries'] == [
                # The annotations have been added to each entry:
                {'pmid': 1, 'foo': 'bar'},
                {'pmid': 2, 'foo': 'baz'}
            ]
        if id_ == 4:
            assert omim_entry['pubmed_entries'] == [
                {'pmid': 1, 'foo': 'bar'},
                {'pmid': 3}  # No annotation added for 3
            ]

    # Test the original dicts were not modified
    for omim_entry in omim_entries:
        for pubmed_entry in omim_entry.get('pubmed_entries', []):
            assert 'foo' not in pubmed_entry

