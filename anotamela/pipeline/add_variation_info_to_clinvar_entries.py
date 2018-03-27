from copy import deepcopy


def add_variation_info_to_clinvar_entries(clinvar_entries, clinvar_variations):
    """
    Given a list of *clinvar_entries* (each one a dictionary), add information
    to each entry from its associated variation (via its "variant_id"). The
    assocaited variation will be looked among the provided *clinvar_variations*.

    Returns a copy of the entries, modified.
    """
    modified_entries = deepcopy(clinvar_entries)

    for clinvar_entry in modified_entries:

        # Find the associated ClinVar Variation:
        variations = [variation for variation in clinvar_variations
                      if variation['variation_id'] == clinvar_entry.get('variant_id')]

        if not variations:
            continue

        variation = variations[0]

        # Add its info to the ClinVar entry
        keys_to_add = [
            'variation_name',
            'variation_type',
        ]
        for key in keys_to_add:
            clinvar_entry[key] = variation[key]

    return modified_entries
