from bs4 import BeautifulSoup

from anotamela.annotators.base_classes import EntrezAnnotator


class DbsnpEntrezAnnotator(EntrezAnnotator):
    """
    Provider Entrez SNP annotations. Usage:

        > dbsnp_entrez_annotator = DbsnpEntrezAnnotator()
        > dbsnp_entrez_annotator.annotate('rs123 rs268'.split())
        # => { 'rs123': ... , 'rs268': ... }

    """
    SOURCE_NAME = 'dbsnp_entrez'
    LINKOUT_NAMES = {'1': 'snp', '5': 'pubmed'}
    ENTREZ_SERVICE = 'epost'
    ENTREZ_PARAMS = {'db': 'snp', 'retmode': 'xml'}

    @staticmethod
    def _annotations_by_id(_, xml_with_many_variants):
        """Splits the XML per variant and returns a dicitonary with the form
        { id-1: xml_fragment-1, id-2: ... }."""
        # The _ ignored argument is the list of IDs, which here is not used
        # because the ID is taken from the xml element itself.
        soup = BeautifulSoup(xml_with_many_variants, 'lxml')
        for xml_element in soup.select('rs'):
            id_ = 'rs' + xml_element['rsid']
            yield id_, str(xml_element)

    @staticmethod
    def _parse_id(id_):
        return id_.replace('rs', '')

    @classmethod
    def _parse_annotation(cls, xml_of_single_variant):
        """Parse the XML of a single rs ID from by Entrez. Return a dict of
        annotations for the given rs."""
        soup = BeautifulSoup(xml_of_single_variant, 'lxml')
        ann = {}

        rs_elements = soup.select('rs')
        assert len(rs_elements) == 1
        e = rs_elements[0]

        ann['rsid'] = 'rs' + e['rsid']
        ann['type'] = e['snpclass']

        seq = [child for child in e.children if child.name == 'sequence']
        assert len(seq) == 1
        seq = seq[0]
        ann['alleles'] = seq.observed.text

        ann['links'] = {}
        for link in e.select('rslinkout'):
            resource = link['resourceid']
            resource_name = cls.LINKOUT_NAMES.get(resource, resource)
            value = link['linkvalue']
            ann['links'][resource_name] = ann['links'].get(resource_name) or []
            ann['links'][resource_name].append(value)

        ann['synonyms'] = []
        for synonym in e.select('mergehistory'):
            ann['synonyms'].append('rs' + synonym['rsid'])

        # I'm removing 'clinical significance' since I found too many
        # inconsistencies between these data and ClinVar's own data on the
        # *same* SNPs. I'd rather not have information than misleading
        # information:

        #  clinsigs = [sig.text
                    #  for sig in e.select('phenotype clinicalsignificance')]
        #  ann['clinical_significance'] = ','.join(clinsigs) or None

        ann['hgvs'] = []
        for hgvs in e.select('hgvs'):
            ann['hgvs'].append(hgvs.text)

        ann['frequency'] = e.frequency and e.frequency.attrs

        ann['fxn'] = [fx.attrs for fx in e.select('fxnset')]

        return ann

