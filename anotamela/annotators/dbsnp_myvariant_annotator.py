from anotamela.annotators.base_classes import MyVariantAnnotator


class DbsnpMyvariantAnnotator(MyVariantAnnotator):
    SOURCE_NAME = 'dbsnp_myvariant'
    SCOPES = 'dbsnp.rsid'
    FIELDS = 'dbsnp'

    @staticmethod
    def _parse_annotation(hits_group):
        dbsnp_data = [hit['dbsnp'] for hit in hits_group if 'dbsnp' in hit]
        if dbsnp_data:
            return dbsnp_data

