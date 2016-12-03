from myvariant import MyVariantInfo

from anotamela.annotators import AnnotatorWithCache, MyVariantAnnotator


class ClinvarRsAnnotator(MyVariantAnnotator, AnnotatorWithCache):
    SOURCE_NAME = 'clinvar_rs'

    def _batch_query(self, ids, _, __):
        mv = MyVariantInfo()
        hits = mv.querymany(ids, scopes='clinvar.rsid', fields='clinvar')
        annotations =  {hit['query']: hit for hit in hits
                        if 'notfound' not in hit}
        self.cache.set(annotations, namespace=self.SOURCE_NAME)
        return annotations

    @staticmethod
    def _parse_annotation(raw_annotation):
        return raw_annotation['clinvar']

