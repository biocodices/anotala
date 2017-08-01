def extract_gwas_traits(gwas_annotations):
    """
    Given a list of gwas_annotations, extract the traits
    of the associations and return them in a list.
    """
    traits = [annotation.get('trait') for annotation in gwas_annotations]
    return sorted({trait for trait in traits if trait})

