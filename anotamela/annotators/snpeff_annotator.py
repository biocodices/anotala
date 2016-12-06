from anotamela.annotators import MyVariantAnnotator


class SnpeffAnnotator(MyVariantAnnotator):
    SOURCE_NAME = 'snpeff_via_myvariant'
    SCOPES = 'dbsnp.rsid'
    FIELDS = 'snpeff'

    @staticmethod
    def _parse_annotation(raw_annotation):
        return raw_annotation['snpeff']

