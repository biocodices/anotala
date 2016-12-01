from os.path import isfile, expanduser
import time
import logging

from Bio import Entrez
from bs4 import BeautifulSoup
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

    You may change the batch size of requests to Entrez by modifying once:

        > DbsnpEntreAnnotator.BATCH_SIZE = 1500  # Default is 1000

    """
    SOURCE_NAME = 'dbsnp_entrez'
    BATCH_SIZE = 1000

    def _batch_query(self, ids, parallel=BATCH_SIZE, sleep_time=0):
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
            soup = BeautifulSoup(xml, 'lxml')
            for rs_element in soup.select('rs'):
                xml_fragment = str(rs_element)
                rs_id = 'rs' + rs_element['rsid']
                group_annotations[rs_id] = xml_fragment

            self.cache.set(group_annotations, namespace=self.SOURCE_NAME)

            annotations.update(group_annotations)
            time.sleep(sleep_time)

        return annotations

    def set_email_for_entrez(self):
        email_filepath = expanduser('~/.mail_address_for_Entrez')

        if not isfile(email_filepath):
            msg = ('Please set a mail for Entrez in {}. Entrez will notify '
                   'that mail before banning you if your usage is too high.')
            raise Exception(msg)

        with open(email_filepath) as f:
            Entrez.email = f.read().strip()
