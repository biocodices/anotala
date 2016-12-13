from anotamela.annotators.base_classes import MyVariantAnnotator


class DbsnpMyvariantAnnotator(MyVariantAnnotator):
    SOURCE_NAME = 'dbsnp_myvariant'
    SCOPES = 'dbsnp.rsid'
    FIELDS = 'dbsnp'

    @staticmethod
    def _parse_hit(hit):
        return hit['dbsnp']

