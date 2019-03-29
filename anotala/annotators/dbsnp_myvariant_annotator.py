from anotala.annotators.base_classes import MyVariantAnnotator


class DbsnpMyvariantAnnotator(MyVariantAnnotator):
    SOURCE_NAME = 'dbsnp_myvariant'
    SCOPES = 'dbsnp.rsid'
    FIELDS = 'dbsnp'

    @staticmethod
    def _parse_hit(hit):
        parsed_hit = hit['dbsnp']
        if 'gene' in parsed_hit and isinstance(parsed_hit['gene'], dict):
            parsed_hit['gene'] = [parsed_hit['gene']]
        return parsed_hit

