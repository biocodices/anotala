from more_itertools import one
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
            variation_id = variation_report['VariationID']
            yield (variation_id, str(variation_report))

    @classmethod
    def _parse_annotation(cls, variation_xml):
        info = {}

        soup = BeautifulSoup(variation_xml, 'lxml-xml')
        variation_reports = soup.select('VariationReport')
        variation_report = one(variation_reports)

        info['variation_id'] = cls._extract_variation_id(variation_report)
        info['variation_name'] = cls._extract_variation_name(variation_report)

        if not info['variation_name'] == 'Multiple Alleles':
            info['variation_type'] = cls._extract_variation_type(variation_report)

        info['genes'] = cls._extract_genes(variation_report)
        info['clinical_assertions'] = \
            cls._extract_clinical_assertions(variation_report)

        return info

    @staticmethod
    def _extract_variation_id(variation_soup):
        """Extract the Variation ID from a ClinVar variation HTML soup."""
        return variation_soup['VariationID']

    @staticmethod
    def _extract_variation_name(variation_soup):
        """Extract the Variation name from a ClinVar variation HTML soup."""
        return variation_soup['VariationName']

    @staticmethod
    def _extract_variation_type(variation_soup):
        """Extract the Variation type from a ClinVar variation HTML soup."""
        return variation_soup['VariationType']

    @staticmethod
    def _extract_clinical_assertions(variation_soup):
        """Extract the clinical assertions from a ClinVar variation."""
        clinical_assertions = []
        selector = 'ClinicalAssertionList > GermlineList > Germline'

        for germline in variation_soup.select(selector):
            info = {}
            info['submitter_name'] = germline['SubmitterName']
            info['date_last_submitted'] = germline['DateLastSubmitted']
            clinsig = one(germline.select('ClinicalSignificance'))
            info['clinical_significance'] = clinsig.select_one('Description').text
            info['method'] = clinsig.select_one('Method').text

            info['phenotypes'] = []
            for pheno in germline.select('PhenotypeList > Phenotype'):
                xrefs = pheno.select_one('XRefList')
                info['phenotypes'].append({
                    'name': pheno['Name'],
                    'omim_id': xrefs.select_one('XRef[DB=OMIM]')['ID'],
                })

            clinical_assertions.append(info)

        return clinical_assertions

    @staticmethod
    def _extract_genes(variation_soup):
        """
        Extract the Variation gene list from a ClinVar variation HTML soup.
        """
        genes = variation_soup.select('GeneList > Gene')

        gene_dicts = []

        for gene in genes:
            gene_info = {
                'symbol': gene.get('Symbol'),
                'full_name': gene.get('FullName'),
                'strand': gene.get('strand'),
                'entrez_id': gene.get('GeneID'),
            }

            if gene.get('HGNCID'):
                gene_info['hgnc_id'] = gene['HGNCID'].replace('HGNC:', '')

            omim = gene.select_one('OMIM')
            if omim:
                gene_info['omim_id'] = omim.text

            gene_dicts.append(gene_info)

        return gene_dicts

