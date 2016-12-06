from anotamela.annotators import MyVariantAnnotator


class DbsnpMyvariantAnnotator(MyVariantAnnotator):
    SOURCE_NAME = 'dbsnp_myvariant'
    SCOPES = 'dbsnp.rsid'
    FIELDS = 'dbsnp'

    @staticmethod
    def _parse_annotation(raw_annotation):
        return raw_annotation['dbsnp']

