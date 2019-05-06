from mygene import MyGeneInfo

from anotala.annotators.base_classes import WebAnnotatorWithCache
from anotala.helpers import grouped


class MygeneAnnotator(WebAnnotatorWithCache):
    """
    Provides gene annotations given one or more gene entrez IDs.
    """
    SOURCE_NAME = 'mygene'
    ANNOTATIONS_ARE_JSON = True
    VERBOSE = False
    TAXID = 9606
    # ^ Human ID, used to avoid annotating with genes from another species
    # It can be replaced at Runtime to annotate with other species
    BATCH_SIZE = 1000

    def _batch_query(self, ids):
        """
        Uses mygene.info service to query many Entrez gene IDs. It returns a
        dict of {id-1: result-1, id-2: ... } with the IDs that were found (i.e.
        leaves out the not found ones).
        """
        if not hasattr(self, 'mg'):
            self.mg = MyGeneInfo()

        for batch_of_ids in grouped(ids, self.BATCH_SIZE):
            batch_annotations = {}
            for hit in self.mg.querymany(batch_of_ids, scopes='entrezgene',
                                         fields='all', verbose=self.VERBOSE):
                if 'notfound' not in hit and hit['taxid'] == self.TAXID:
                    batch_annotations[hit['query']] = hit
            yield batch_annotations

    @staticmethod
    def _parse_annotation(raw_annotation):
        fields_to_keep = ('name symbol entrezgene MIM HGNC uniprot '
                          'type_of_gene').split()
        annotation = {k: v for k, v in raw_annotation.items()
                      if k in fields_to_keep}
        if 'uniprot' in annotation:
            swissprot_id = annotation['uniprot'].get('Swiss-Prot')
            if swissprot_id:
                annotation['swissprot'] = swissprot_id
            del(annotation['uniprot'])
        return annotation

