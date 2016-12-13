from anotamela.annotators.base_classes import MyVariantAnnotator


class SnpeffAnnotator(MyVariantAnnotator):
    SOURCE_NAME = 'snpeff_myvariant'
    SCOPES = 'dbsnp.rsid'
    FIELDS = 'snpeff'

    @staticmethod
    def _parse_hit(hit):
        return hit['snpeff']

