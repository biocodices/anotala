from anotamela.pipeline import add_variation_info_to_clinvar_entries


def test_add_variation_info_to_clinvar_entries():
    clinvar_entries = [
        {'variant_id': '1'},
        {'variant_id': None},
        {}
    ]
    clinvar_variations = [
        {
            'variation_id': '1',
            'variation_name': 'Name-1',
            'variation_type': 'Type-1',
        },
        {
            'variation_id': '2',
            'variation_name': 'Name-2',
            'variation_type': 'Type-2',
        },
    ]

    modified_entries = add_variation_info_to_clinvar_entries(clinvar_entries,
                                                             clinvar_variations)
    assert modified_entries[0]['variation_name'] == 'Name-1'
    assert modified_entries[0]['variation_type'] == 'Type-1'
