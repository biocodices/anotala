from itertools import groupby
from operator import itemgetter

from mygene import MyGeneInfo
import numpy as np

from anotamela.annotators import AnnotatorWithCache


class MygeneAnnotator(AnnotatorWithCache):
    """
    Provides gene annotations given one or more gene symbols (e.g. TTN).
    """
    SOURCE_NAME = 'mygene'
    ANNOTATIONS_ARE_JSON = True

    def annotate(self, ids, use_cache=True, use_web=True, parse_data=True):
        # This wrapper around AnnotatorWithCache.annotate is meant to remove
        # the parallel and sleep_time arguments. Those options are superflous
        # since the mygene client already deals with batch queries.
        return super().annotate(ids, parallel=len(ids), sleep_time=0,
                                use_cache=use_cache, use_web=use_web,
                                parse_data=parse_data)

    def _batch_query_and_cache(self, ids, _, __):
        """
        Uses mygene.info service to query many entrez gene IDs. It returns a
        dict of {id: result, ... } with the IDs that were found (i.e. leaves
        out the not found ones). The successful results are also cached.

        The ignored arguments _ and __ are there to handle the <parallel>
        and <sleep_time> arguments that AnnotatorWithCache.annotate
        assumes, since mygene client already handles batch queries.
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
        del(annotation['uniprot']['TrEMBL'])
        return annotation

