from collections import defaultdict


def group_swissprot_variants_by_rsid(swissprot_variants):
    """
    Given a list of SwissProt variants, group them in a dict with the rsids
    as keys and the corresponding list of entries for each rsid as values.
    """
    rs_to_swissprot_variants = defaultdict(list)

    for swissprot_variant in swissprot_variants:
        rsid = swissprot_variant.get('rsid')
        if rsid:
            rs_to_swissprot_variants[rsid].append(swissprot_variant)

    return dict(rs_to_swissprot_variants)

