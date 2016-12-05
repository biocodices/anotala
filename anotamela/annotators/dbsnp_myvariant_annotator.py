from anotamela.annotators import AnnotatorWithCache, MyVariantAnnotator


class DbsnpMyvariantAnnotator(MyVariantAnnotator, AnnotatorWithCache):
    SOURCE_NAME = 'dbsnp_myvariant'
    SCOPES = 'dbsnp.rsid'
    FIELDS = 'dbsnp'

    @staticmethod
    def _parse_annotation(raw_annotation):
        return raw_annotation['dbsnp']

