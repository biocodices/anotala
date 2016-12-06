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
    BATCH_SIZE = 1000

    def _batch_query_and_cache(self, ids, _, __):
        # The ignored arguments _ and _ _ are there to handle the <parallel>
        # and <sleep_time> arguments that AnnotatorWithCache.annotate assumes.
        # The Entrez ePOST service already handles batch queries.
        ids = [id_.replace('rs', '') for id_ in ids]
        entrez_helper = EntrezHelper()
        entrez_params = {
                'db_name': 'snp',
                'ids': ids,
                'rettype': 'xml',
                'batch_size': self.BATCH_SIZE,
                'xml_element_tag': 'rs',
                'xml_id_attribute': 'rsid'
            }
        annotations = {}
        for batch_annotations in entrez_helper.post_query(**entrez_params):
            annotations.update(batch_annotations)
            self.cache.set(batch_annotations, namespace=self.SOURCE_NAME)
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

