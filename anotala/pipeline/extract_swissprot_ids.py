def extract_swissprot_ids(mygene_annotations):
    """
    Given a pandas Series of mygene annotations, return a list of the unique
    Swissprot IDs in them.
    """
    uniprot_ids = set()

    for gene in mygene_annotations.fillna(False):
        if gene and 'swissprot' in gene:
            ids = gene['swissprot']
            if not isinstance(ids, list):
                ids = [ids]
            uniprot_ids = uniprot_ids | set(ids)

    return list(uniprot_ids)

