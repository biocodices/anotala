from mygene import MyGeneInfo

from anotamela.annotators.base_classes import AnnotatorWithCache


class MygeneAnnotator(AnnotatorWithCache):
    """
    Provides gene annotations given one or more gene symbols (e.g. TTN).
    """
    SOURCE_NAME = 'mygene'
    ANNOTATIONS_ARE_JSON = True
    TAXID = 9606  # Human taxon ID, used to annotate the gene of this species

    def _batch_query(self, ids):
        """
        Uses mygene.info service to query many Entrez gene IDs. It returns a
        dict of {id-1: result-1, id-2: ... } with the IDs that were found (i.e.
        leaves out the not found ones).
        """
        if not hasattr(self, 'mg'):
            self.mg = MyGeneInfo()

        for hit in self.mg.querymany(ids, fields='all'):
            if 'notfound' not in hit and hit['taxid'] == self.TAXID:
                yield hit['query'], hit

    @staticmethod
    def _parse_annotation(raw_annotation):
        fields_to_keep = ('name symbol entrezgene MIM HGNC uniprot '
                          'type_of_gene').split()
        annotation = {k: v for k, v in raw_annotation.items()
                      if k in fields_to_keep}
        annotation['swissprot'] = annotation['uniprot']['Swiss-Prot']
        del(annotation['uniprot'])
        return annotation

