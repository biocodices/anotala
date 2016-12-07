from mygene import MyGeneInfo

from anotamela.annotators.base_classes import AnnotatorWithCache


class MygeneAnnotator(AnnotatorWithCache):
    """
    Provides gene annotations given one or more gene symbols (e.g. TTN).
    """
    SOURCE_NAME = 'mygene'
    ANNOTATIONS_ARE_JSON = True

    def _batch_query_and_cache(self, ids):
        """
        Uses mygene.info service to query many entrez gene IDs. It returns a
        dict of {id: result, ... } with the IDs that were found (i.e. leaves
        out the not found ones). The successful results are also cached.
        """
        if not hasattr(self, 'mg'):
            self.mg = MyGeneInfo()
        hits = self.mg.querymany(ids, fields='all')
        human_taxid = 9606
        annotations = {hit['query']: hit for hit in hits
                       if 'notfound' not in hit and hit['taxid'] == human_taxid}
        self.cache.set(annotations, namespace=self.SOURCE_NAME)
        return annotations

    @staticmethod
    def _parse_annotation(raw_annotation):
        fields_to_keep = ('name symbol entrezgene MIM HGNC uniprot '
                          'type_of_gene').split()
        annotation = {k: v for k, v in raw_annotation.items()
                      if k in fields_to_keep}
        annotation['swissprot'] = annotation['uniprot']['Swiss-Prot']
        del(annotation['uniprot'])
        return annotation

