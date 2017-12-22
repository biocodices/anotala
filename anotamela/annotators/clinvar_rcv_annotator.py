from more_itertools import one
from bs4 import BeautifulSoup

from anotamela.annotators.base_classes import EntrezAnnotator



class ClinvarRCVAnnotator(EntrezAnnotator):
    """
    Annotates ClinVar RCV accessions, using Biopython's Entrez service.
    """
    SOURCE_NAME = 'clinvar'
    ENTREZ_SERVICE = 'efetch'
    ENTREZ_PARAMS = {
        'db': 'clinvar',
        'rettype': 'clinvarset',
    }

    @staticmethod
    def _annotations_by_id(ids, multi_accession_xml):
        """
        Given an XML response with many <ClinVarSet> elements,
        yield tuples of (RCV accession number, ClinVarSet XML).
        """
        soup = BeautifulSoup(multi_accession_xml, 'lxml-xml')

        for clinvar_set in soup.select('ClinVarSet'):
            accessions = [element['Acc'] for element in
                          clinvar_set.find_all('ClinVarAccession')
                          if element.get('Type') == 'RCV']

            try:
                accession = one(accessions)
            except ValueError:
                raise ValueError('{} accessions found (expecting one).'
                                 .format(len(accessions)))

            yield (accession, str(clinvar_set))

    @classmethod
    def _parse_annotation(cls, clinvar_set_xml):
        info = {}

        soup = BeautifulSoup(clinvar_set_xml, 'lxml-xml')

        info['accession'] = cls._extract_accession(soup)
        info['entry_type'] = cls._extract_entry_type(soup)
        info['title'] = cls._extract_title(soup)
        info['attributes'] = cls._extract_attributes(soup)

        return info

    @staticmethod
    def _extract_title(soup):
        return soup.find('Title').text.strip()

    @staticmethod
    def _extract_attributes(soup):
        attributes = []
        for attribute in soup.select('MeasureSet AttributeSet Attribute'):
            info = {key.lower(): val for key, val in attribute.attrs.items()}
            info['full_name'] = attribute.text.strip()
            attributes.append(info)
        return attributes

    @staticmethod
    def _extract_dbsnp_id(soup):
        xref = soup.find('XRef', attrs={'DB': 'dbSNP', 'Type': 'rs'})
        return 'rs' + xref['ID']

    @staticmethod
    def _extract_accession(soup):
        return soup.find('ClinVarAccession', attrs={'Type': 'RCV'})['Acc']

    @staticmethod
    def _extract_entry_type(soup):
        return soup.find('MeasureSet')['Type']

