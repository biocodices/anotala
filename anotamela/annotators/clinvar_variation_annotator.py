from collections import defaultdict, Counter

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
        clinical_assertions = cls._extract_clinical_assertions(variation_report)
        info['clinical_assertions'] = clinical_assertions
        info['clinical_summary'] = \
            dict(cls._generate_clinical_summary(clinical_assertions))
        info['clinical_significance'] = \
            cls._extract_clinical_significance(variation_report)
        info['associated_phenotypes'] = \
            cls._associated_phenotypes(clinical_assertions)

        info['alleles'] = cls._extract_alleles(variation_report)

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
        return variation_soup.get('VariationType')

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
                pheno_info = {'name': pheno['Name']}
                omim = xrefs.select_one('XRef[DB=OMIM]')
                if omim:
                    pheno_info['omim_id'] = omim['ID']

                info['phenotypes'].append(pheno_info)


            clinical_assertions.append(info)

        return clinical_assertions

    @staticmethod
    def _extract_clinical_significance(variation_soup):
        """Extract the single "Clinical significance" that describes a
        CliVar variation. It might be labeled as "Conflicting" if different
        clinical assertions do not agree."""
        selector = 'ObservationList > ClinicalSignificance > Description'
        clin_sig = variation_soup.select_one(selector)
        if clin_sig:
            return clin_sig.text


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
            info.update(cls._extract_xrefs(allele))
            info['consequences'] = cls._extract_molecular_consequences(allele)
            info['frequencies'] = cls._extract_allele_frequencies(allele)

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

        if g37:
            # Copy Numbers have an "innerStart" instead of "start"
            start = g37.get('start') or g37.get('innerStart')
            info['start_g37'] = int(start)

            # Copy Numbers have an "innerStop" instead of "stop"
            stop = g37.get('stop') or g37.get('innerStop')
            info['stop_g37'] = int(stop)

            info['accession_g37'] = g37.get('Accession')

            length = g37.get('variantLength')
            if length:
                info['length_g37'] = int(length)

            info['ref_g37'] = g37.get('referenceAllele')
            info['alt_g37'] = g37.get('alternateAllele')
            info['chrom_g37'] = g37.get('Chr')

        if g38:
            start = g38.get('start') or g38.get('innerStart')
            info['start_g38'] = int(start)

            stop = g38.get('stop') or g38.get('innerStop')
            info['stop_g38'] = int(stop)

            info['accession_g38'] = g38.get('Accession')

            length = g38.get('variantLength')
            if length:
                info['length_g38'] = int(length)

            info['ref_g38'] = g38.get('referenceAllele')
            info['alt_g38'] = g38.get('alternateAllele')
            info['chrom_g38'] = g38.get('Chr')

        return info

    @staticmethod
    def _extract_allele_hgvs(allele):
        """Given an <Allele> BS node, extract genomic, coding, and protein
        HGVS changes."""
        info = {}
        hgvs = allele.select('HGVSlist')

        if hgvs:
            hgvs = one(hgvs)
        else:
            return info

        g37 = hgvs.select_one('HGVS[Assembly=GRCh37]')
        if g37:
            info['genomic_change_g37'] = g37.get('Change')
            info['genomic_change_g37_accession'] = g37.get('AccessionVersion')
            info['genomic_change_g37_name'] = g37.text

        g38 = hgvs.select_one('HGVS[Assembly=GRCh38]')
        if g38:
            info['genomic_change_g38'] = g38.get('Change')
            info['genomic_change_g38_accession'] = g38.get('AccessionVersion')
            info['genomic_change_g38_name'] = g38.text

        c = hgvs.find('HGVS', attrs={'Type': 'HGVS, coding, RefSeq'})
        if c:
            info['coding_change'] = c.get('Change')
            info['coding_change_accession'] = c.get('AccessionVersion')

        p = hgvs.find('HGVS', attrs={'Type': 'HGVS, protein, RefSeq'})
        if p:
            info['protein_change'] = p.get('Change')
            info['protein_change_accession'] = p.get('AccessionVersion')

        return info

    @staticmethod
    def _extract_xrefs(allele):
        """Given an <Allele> BS node, extract the external DB references."""
        info = {}

        dbsnp = allele.find('XRef', attrs={'DB': 'dbSNP', 'Type': 'rs'})
        if dbsnp:
            info['dbsnp_id'] = '{}{}'.format(dbsnp['Type'], dbsnp['ID'])

        omim = allele.find('XRef', attrs={'DB': 'OMIM'})
        if omim:
            info['omim_id'] = omim['ID']

        uniprot = allele.find('XRef', attrs={'DB': 'UniProtKB'})
        if uniprot:
            info['uniprot_id'] = uniprot['ID']

        return info

    @staticmethod
    def _extract_molecular_consequences(allele):
        """Given an <Allele> BS node, extract the molecular consequences."""
        consequences = []
        for consequence in allele.select('MolecularConsequence'):
            info = {
                'hgvs': consequence['HGVS'],
                'function': consequence['Function'],
            }
            consequences.append(info)

        return consequences

    @staticmethod
    def _extract_allele_frequencies(allele):
        freq_per_allele = defaultdict(dict)

        for frequency in allele.select('AlleleFrequency'):
            allele = frequency['MinorAllele']
            source = frequency['Type']
            value = float(frequency['Value'])

            freq_per_allele[allele][source] = value

        return dict(freq_per_allele)

    @staticmethod
    def _generate_clinical_summary(clinical_assertions):
        return Counter(assertion['clinical_significance']
                       for assertion in clinical_assertions)

    @staticmethod
    def _associated_phenotypes(clinical_assertions):
        phenotype_names = [phenotype['name']
                           for clinical_assertion in clinical_assertions
                           for phenotype in clinical_assertion['phenotypes']]
        return sorted(set(phenotype_names))

