from anotamela.annotators.base_classes import ParallelAnnotator


class GwasCatalogAnnotator(ParallelAnnotator):
    SOURCE_NAME = 'gwas_catalog'
    ANNOTATIONS_ARE_JSON = True

    @staticmethod
    def _url(id_):
        # Gwas Catalog will take any search term here, but we typically
        # use it for rs IDs:
        url = 'https://www.ebi.ac.uk/gwas/api/search?q="{}"'.format(id_)

        # This grouping is optional, but facilitates the parsing afterwards:
        url += '&group=true&group.by=resourcename'
        return url

    @staticmethod
    def _parse_annotation(response):
        # This parsing depends on the URL above having the grouping options:
        groups = {group['groupValue']: group['doclist']['docs']
                  for group in response['grouped']['resourcename']['groups']}

        if not groups:
            return

        assert 'study' in groups
        assert 'diseasetrait' in groups

        return groups['association']

