from anotamela.annotators.base_classes import MyVariantAnnotator


class DbnsfpAnnotator(MyVariantAnnotator):
    SOURCE_NAME = 'dbnsfp'
    SCOPES = 'dbsnp.rsid'
    FIELDS = 'dbnsfp'

    @staticmethod
    def _parse_hit(hit):
        return hit.get('dbnsfp')

