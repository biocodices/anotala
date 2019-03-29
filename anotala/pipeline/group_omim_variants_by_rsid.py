from collections import defaultdict


def group_omim_variants_by_rsid(omim_variants):
    """
    Given a list of omim variant entries, group them in a dict with the rsids
    as keys and the corresponding list of entries for each rsid as values.
    """
    # There might be more than one OMIM variant for a given rs --the typical
    # case is the one of multiallelic SNPs. Hence, we return the OMIM
    # annotation for an rs always as a *list* of entries:
    rs_to_omim_variants = defaultdict(list)

    for omim_variant in omim_variants:
        # One OMIM variant might also correspond to many rs IDs, so we assign
        # the entry to all of them. That means one entry might be repeated
        # across different rsids:
        for rs in omim_variant.get('rsids', []):
            rs_to_omim_variants[rs].append(omim_variant)

    return dict(rs_to_omim_variants)

