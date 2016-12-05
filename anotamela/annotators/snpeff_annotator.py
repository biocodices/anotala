from anotamela.annotators import AnnotatorWithCache, MyVariantAnnotator


class SnpeffAnnotator(MyVariantAnnotator, AnnotatorWithCache):
    SOURCE_NAME = 'snpeff_via_myvariant'
    SCOPES = 'dbsnp.rsid'
    FIELDS = 'snpeff'

    @staticmethod
    def _parse_annotation(raw_annotation):
        return raw_annotation['snpeff']

