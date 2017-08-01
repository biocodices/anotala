def extract_ensembl_consequence(ensembl_annotation):
    """
    Extract the 'most_severe_consequence' from a dictionary of ensembl
    annotations for a given variant.
    """
    if not ensembl_annotation:
        return None

    consequence = ensembl_annotation.get('most_severe_consequence')

    if not consequence or str(consequence) == 'nan':
        return None

    return consequence

