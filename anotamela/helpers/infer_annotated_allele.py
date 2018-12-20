from anotamela.helpers.variation_patterns import (
    SNP_RE,
    SYN_SNP_RE,
    INDEL_RE,
)


def infer_annotated_allele(mutation):
    """Given a *mutation* like 'c.123A>G', infer the allele 'G'."""
    if not mutation:
        return

    allele = mutation  # If all fails, return the same mutation

    for regex in [SNP_RE, SYN_SNP_RE, INDEL_RE]:
        match = regex.search(mutation)
        if match:
            allele = match.group('new_allele')

    return allele
