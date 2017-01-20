def extract_entrez_genes(dbsnp_entries, field='geneid'):
    """
    Given a list of Myvariant *dbsnp_entries* as annotated by
    DbsnpMyvariantAnnotator, extract the Entrez Gene 'geneid' or 'symbol'.

    Returns a list of geneids or symbols, according to the passed *field*.
    """
    if not dbsnp_entries:
        return []

    values = {gene.get(field) for entry in dbsnp_entries
                              for gene in entry.get('gene', {})}
    return list(values)

