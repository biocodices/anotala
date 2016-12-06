from anotamela.annotators import AnnotatorWithCache, EntrezAnnotator


class GeneEntrezAnnotator(EntrezAnnotator, AnnotatorWithCache):
    """
    Provider of Entrez Gene annotations. XML responses are cached and then
    parsed. Usage:

        > gene_entrez = GeneEntrezAnnotator()
        > gene_entrez.annotate('1756')  # Also accepts a list of IDs
        # => { '1756': ... }
    """
    SOURCE_NAME = 'gene_entrez'
    ANNOTATIONS_ARE_JSON = True
    ENTREZ_PARAMS = {'db': 'gene'}
    USE_ENTREZ_READER = True

    @staticmethod
    def _annotations_by_id(multigene_response):
        """For each element in the list that Entrez returns, find which is the
        gene ID for that element and return a tuple (gene_id, element)."""
        for gene in multigene_response:
            try:
                locus_id = [key for key in gene['Entrezgene_unique-keys']
                            if key['Dbtag_db'] == 'LocusID']
                assert len(locus_id) == 1
                gene_id = locus_id[0]['Dbtag_tag']['Object-id']['Object-id_id']
                yield (gene_id, gene)
            except Exception:
                import ipdb; ipdb.set_trace()
                raise

