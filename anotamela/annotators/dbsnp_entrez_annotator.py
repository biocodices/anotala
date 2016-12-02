from os.path import isfile, expanduser
import time
import logging

from Bio import Entrez
from tqdm import tqdm

from anotamela.annotators import AnnotatorWithCache


logger = logging.getLogger(__name__)


class DbsnpEntrezAnnotator(AnnotatorWithCache):
    """
    Provider of DbSNP annotations taken from the Entrez service. Responses are
    cached. Usage:

        > dbsnp_entrez_annotator = DbsnpEntrezAnnotator()
        > dbsnp_entrez_annotator.annotate('rs123 rs268'.split())
        # => { 'rs123': ... , 'rs268': ... }
    """
    SOURCE_NAME = 'dbsnp_entrez'
    ANNOTATION_TYPE = 'XML'
    LINKOUT_NAMES = {'1': 'snp', '5': 'pubmed'}

    def _batch_query(self, ids, parallel, sleep_time):
        if not Entrez.email:
            self.set_email_for_entrez

        # Entrez POST queries are a two step process. You first POST your query
        # and get a WebEnv identifier and a QueryKey. Then you query using
        # those IDs in batches. More info:
        # http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec:entrez-webenv
        logger.info('Create a serverside Job and get its ID')
        ids = [id_.replace('rs', '') for id_ in ids]
        handle = Entrez.epost(db='snp', id=','.join(ids))
        job = Entrez.read(handle)

        total = len(ids)
        msg = 'Fetch {} entries from Entrez DbSNP in batches of {}'
        logger.info(msg.format(total, parallel))

        annotations = {}
        for offset in tqdm(list(range(0, total, parallel))):
            # Fetch the XML
            fetch_handle = Entrez.efetch(db='snp', rettype='xml',
                    retstart=offset, retmax=parallel,
                    webenv=job['WebEnv'], query_key=job['QueryKey'])
            xml = fetch_handle.read()
            fetch_handle.close()

            # Split the XML per variant
            group_annotations = {}
            soup = self._make_xml_soup(xml)
            for rs_element in soup.select('rs'):
                xml_fragment = str(rs_element)
                rs_id = 'rs' + rs_element['rsid']
                group_annotations[rs_id] = xml_fragment

            self.cache.set(group_annotations, namespace=self.SOURCE_NAME)

            annotations.update(group_annotations)
            time.sleep(sleep_time)

        return annotations

    @classmethod
    def _parse_annotation(cls, annotation):
        """
        Parse an XML for a single rs ID from by Entrez. Return a dict.
        """
        ann = {}
        soup = cls._make_xml_soup(annotation)

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

        ann['frequency'] = e.frequency.attrs
        ann['fxn'] = [fx.attrs for fx in e.select('fxnset')]

        return ann

    def set_email_for_entrez(self):
        email_filepath = expanduser('~/.mail_address_for_Entrez')

        if not isfile(email_filepath):
            msg = ('Please set a mail for Entrez in {}. Entrez will notify '
                   'that mail before banning you if your usage is too high.')
            raise Exception(msg)

        with open(email_filepath) as f:
            Entrez.email = f.read().strip()
