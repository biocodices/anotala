import logging

from anotamela.annotators import AnnotatorWithCache
from anotamela.helpers import make_xml_soup, EntrezHelper


logger = logging.getLogger(__name__)


class DbsnpEntrezAnnotator(AnnotatorWithCache):
    """
    Provider of DbSNP annotations taken from the Entrez service. XML responses
    are cached and then parsed. Usage:

        > dbsnp_entrez_annotator = DbsnpEntrezAnnotator()
        > dbsnp_entrez_annotator.annotate('rs123 rs268'.split())
        # => { 'rs123': ... , 'rs268': ... }
    """
    SOURCE_NAME = 'dbsnp_entrez'
    LINKOUT_NAMES = {'1': 'snp', '5': 'pubmed'}

    def _batch_query(self, ids, parallel, _):
        ids = [id_.replace('rs', '') for id_ in ids]
        entrez_params = {
                'db_name': 'snp',
                'ids': ids,
                'rettype': 'xml',
                'batch_size': parallel,
                'xml_element_tag': 'rs',
                'xml_id_attribute': 'rsid'
            }
        entrez_helper = EntrezHelper()
        results = entrez_helper.post_query(**entrez_params)
        annotations = {}
        for id_, annotation in results:
            self.cache.set({id_: annotation}, namespace=self.SOURCE_NAME)
            annotations[id_] = annotation
        return annotations

    @classmethod
    def _parse_annotation(cls, annotation):
        """Parse the XML of a single rs ID from by Entrez. Return a dict of
        annotations for the given rs."""
        ann = {}
        soup = make_xml_soup(annotation)

        rs_elements = soup.select('rs')
        assert len(rs_elements) == 1
        e = rs_elements[0]

        ann['rs_id'] = 'rs' + e['rsid']
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

        clinsigs = [sig.text
                    for sig in e.select('phenotype clinicalsignificance')]
        ann['clinical_significance'] = ','.join(clinsigs) or None

        ann['hgvs'] = []
        for hgvs in e.select('hgvs'):
            ann['hgvs'].append(hgvs.text)

        ann['frequency'] = e.frequency and e.frequency.attrs
        ann['fxn'] = [fx.attrs for fx in e.select('fxnset')]

        return ann

