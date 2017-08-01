def extract_ensembl_consequence(ensembl_annotation):
    if not ensembl_annotation:
        return None

    consequence = ensembl_annotation.get('most_severe_consequence')

    if not consequence or str(consequence) == 'nan':
        return None

    return consequence

