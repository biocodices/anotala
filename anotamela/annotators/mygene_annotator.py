from itertools import groupby
from operator import itemgetter

from mygene import MyGeneInfo
import numpy as np

from anotamela.annotators import AnnotatorWithCache


class MygeneAnnotator(AnnotatorWithCache):
    """
    Provides gene annotations given one or more gene symbols (e.g. TTN).
    """
    ANNOTATIONS_ARE_JSON = True

    def annotate(self, ids, use_cache=True, use_web=True, parse_data=True):
        # This wrapper around AnnotatorWithCache.annotate is meant to remove
        # the parallel and sleep_time arguments. Those options are superflous
        # since the mygene client already deals with batch queries.
        return super().annotate(ids, parallel=len(ids), sleep_time=0,
                                use_cache=use_cache, use_web=use_web,
                                parse_data=parse_data)

    def _batch_query_and_cache(self, ids, _, __):
        # The ignored arguments _ and __ are there to handle the <parallel>
        # and <sleep_time> arguments that AnnotatorWithCache.annotate
        # assumes.
        # Uses mygene.info service to query many IDs across the given
        # scopes. It returns a dict of {id: result, ... } with the IDs that
        # were found (i.e. leaves out the not found ones). The successful
        # results are also cached.
        if not hasattr(self, 'mg'):
            self.mg = MyGeneInfo()
        hits = self.mg.querymany(ids, scopes='symbol', fields='all')
        sort_func = itemgetter('query')
        hits_by_gene = groupby(sorted(results, key=sort_func), sort_func)


        ### SEGUIR ACÁ! max(results, func=score.....)


        annotations = {hit['query']: hit for hit in hits
                       if 'notfound' not in hit}
        # Deal with the multiple hits per gene symbol
        self.cache.set(annotations, namespace=self.SOURCE_NAME)
        return annotations

    @statimethod
    def _parse_annotation(raw_annotation):
        fields_to_keep = ('name summary symbol entrezgene MIM HGNC '
                          'type_of_gene uniprot').split()
        return {k: v for k, v in raw_annotation.items() if k in fields_to_keep}

