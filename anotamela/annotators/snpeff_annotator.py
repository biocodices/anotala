from anotamela.annotators.base_classes import MyVariantAnnotator


class SnpeffAnnotator(MyVariantAnnotator):
    SOURCE_NAME = 'snpeff_myvariant'
    SCOPES = 'dbsnp.rsid'
    FIELDS = 'snpeff'

    @staticmethod
    def _parse_hit(hit):
        entry = hit['snpeff']

        # Make sure the annotations are always a *list* of dictionaries
        # Right now, myvariant client sometimes returns a list and sometimes
        # a single dictionary.
        if not isinstance(entry['ann'], list):
            entry['ann'] = [entry['ann']]

        return entry

