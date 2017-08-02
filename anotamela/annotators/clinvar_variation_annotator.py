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

    @classmethod
    def _extract_alleles(cls, variation_soup):
        """Extract the Variation alleles from a ClinVar variation."""
        alleles = []
        for allele in variation_soup.select('Allele'):
            info = {}
            info.update(cls._extract_allele_basic_info(allele))
            info.update(cls._extract_sequence_info_from_allele(allele))
            info.update(cls._extract_allele_hgvs(allele))

            alleles.append(info)

        return alleles


    @staticmethod
    def _extract_allele_basic_info(allele):
        """Given an <Allele> BS node, extract its basic info into a dict."""
        info = {}
        info['allele_id'] = allele['AlleleID']
        info['name'] = allele.select_one('Name').text
        info['variant_type'] = allele.select_one('VariantType').text
        return info


    @staticmethod
    def _extract_sequence_info_from_allele(allele):
        """Given an <Allele> BS node, extract SequenceLocation info into a dict."""
        g37 = allele.select_one('SequenceLocation[Assembly="GRCh37"]')
        g38 = allele.select_one('SequenceLocation[Assembly="GRCh38"]')

        info = {}

        info['start_g37'] = int(g37['start'])
        info['stop_g37'] = int(g37['stop'])
        info['accession_g37'] = g37['Accession']
        info['length_g37'] = int(g37['variantLength'])
        info['ref_g37'] = g37['referenceAllele']
        info['alt_g37'] = g37['alternateAllele']
        info['chrom_g37'] = g37['Chr']

        info['start_g38'] = int(g38['start'])
        info['stop_g38'] = int(g38['stop'])
        info['accession_g38'] = g38['Accession']
        info['length_g38'] = int(g38['variantLength'])
        info['ref_g38'] = g38['referenceAllele']
        info['alt_g38'] = g38['alternateAllele']
        info['chrom_g38'] = g38['Chr']

        return info


    @staticmethod
    def _extract_allele_hgvs(allele):
        """Given an <Allele> BS node, extract genomic, coding, and protein
        HGVS changes."""
        info = {}
        hgvs = one(allele.select('HGVSlist'))

        g37 = hgvs.select_one('HGVS[Assembly=GRCh37]')
        g38 = hgvs.select_one('HGVS[Assembly=GRCh38]')
        # Can't use CSS selectors when the attribute has whitespace in it.
        c = hgvs.find('HGVS', attrs={'Type': 'HGVS, coding, RefSeq'})
        p = hgvs.find('HGVS', attrs={'Type': 'HGVS, protein, RefSeq'})

        info['genomic_change_g37'] = g37['Change']
        info['genomic_change_g37_accession'] = g37['AccessionVersion']
        info['genomic_change_g38'] = g38['Change']
        info['genomic_change_g38_accession'] = g38['AccessionVersion']
        info['coding_change'] = c['Change']
        info['coding_change_accession'] = c['AccessionVersion']
        info['protein_change'] = p['Change']
        info['protein_change_accession'] = p['AccessionVersion']

        return info
