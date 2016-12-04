from anotamela.annotators import AnnotatorWithCache, MyVariantAnnotator


class ClinvarRsAnnotator(MyVariantAnnotator, AnnotatorWithCache):
    SOURCE_NAME = 'clinvar_rs'

    def _batch_query(self, ids, _, __):
        return self._myvariant_query_and_cache(ids, scopes='clinvar.rsid',
                                               fields='clinvar')

    @staticmethod
    def _parse_annotation(raw_annotation):
        return raw_annotation['clinvar']

