from bs4 import BeautifulSoup

from anotamela.annotators.base_classes import EntrezAnnotator


class ClinvarVariationAnnotator(EntrezAnnotator):
    """
    Annotates Clinvar Variation IDs using Biopython's Entrez service.
    """
    SOURCE_NAME = 'clinvar'
    ENTREZ_SERVICE = 'efetch'
    ENTREZ_PARAMS = {
        'db': 'clinvar',
        'rettype': 'variation',
    }

    @staticmethod
    def _annotations_by_id(ids, xml):
        """
        Given an XML response with many <VariationReport> elements,
        yield tuples of (Variation ID, XML for that variation).
        """
        soup = BeautifulSoup(xml, 'lxml-xml')

        for variation_report in soup.select('VariationReport'):
            variation_id = variation_report['VariationId']
            yield (variation_id, str(variation_report))

    @classmethod
    def _parse_annotation(cls, variation_xml):
        info = {}

        soup = BeautifulSoup(variation_xml, 'lxml-xml')

        info['variation_id'] = cls._extract_variation_id(soup)

        return info

    @staticmethod
    def _extract_variation_id(variation_soup):
        """Extract the Variation ID from a ClinVar variation HTML soup."""
        return variation_soup['VariationId']

    @staticmethod
    def _extract_variation_name(variation_soup):
        """Extract the Variation name from a ClinVar variation HTML soup."""
        return variation_soup['VariationName']

    @staticmethod
    def _extract_variation_type(variation_soup):
        """Extract the Variation type from a ClinVar variation HTML soup."""
        return variation_soup['VariationType']

    @staticmethod
    def _extract_genes(variation_soup):
        """
        Extract the Variation gene list from a ClinVar variation HTML soup.
        """
        genes = variation_soup.select_one('GeneList > Gene')
        # return


