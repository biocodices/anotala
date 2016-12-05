from anotamela.annotators import AnnotatorWithCache, MyVariantAnnotator


class ClinvarRsAnnotator(MyVariantAnnotator, AnnotatorWithCache):
    SOURCE_NAME = 'clinvar_rs'
    SCOPES = 'clinvar.rsid'
    FIELDS = 'clinvar'

    @staticmethod
    def _parse_annotation(raw_annotation):
        return raw_annotation['clinvar']

