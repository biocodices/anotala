import time
from os.path import isfile, expanduser
import logging

from Bio import Entrez
from tqdm import tqdm

from anotamela.helpers import make_xml_soup


logger = logging.getLogger(__name__)


class EntrezHelper():
    def post_query(self, ids, db_name, batch_size, rettype='json',
                   xml_element_tag=None, xml_id_attribute=None):
        """
        Use Entrez POST query service to fetch a list of IDs in the given
        database. Yields tuples of (id, annotation).

        If rettype='xml', set xml_element_tag and xml_id_attribute to let the
        XML parser know which XML elements are associated to each ID. For
        instance, for the db 'snp', set xml_element_tag='rs' and
        xml_id_attribute='rsid'. You have to manually explore the XML once
        to know theese values.
        """
        if not Entrez.email:
            self.set_email_for_entrez()

        # Entrez POST queries are a two step process. You first POST your query
        # and get a WebEnv identifier and a QueryKey.
        logger.info('Create a serverside Job and get its ID')
        handle = Entrez.epost(db=db_name, id=','.join(ids))
        job = Entrez.read(handle)

        # Then you do a second query using those IDs, and you get the results
        # in batches. More info:
        # http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec:entrez-webenv
        total = len(ids)
        msg = 'Fetch {} entries from Entrez DbSNP in batches of {}'
        logger.info(msg.format(total, batch_size))
        for offset in tqdm(list(range(0, total, batch_size))):
            # Fetch the raw response (XML or JSON)
            fetch_handle = Entrez.efetch(db=db_name, rettype=rettype,
                    retstart=offset, retmax=batch_size,
                    webenv=job['WebEnv'], query_key=job['QueryKey'])
            raw_response = fetch_handle.read()
            fetch_handle.close()

            if rettype == 'xml':
                batch_annotations = self.extract_id_annotation_from_xml(
                    raw_response, xml_element_tag, xml_id_attribute, db_name)
                yield dict(batch_annotations)
            else:
                msg = 'No extraction method for rettype={}'.format(rettype)
                raise NotImplementedError(msg)

    def extract_id_annotation_from_xml(self, raw_response, xml_element_tag,
                                       xml_id_attribute, db_name):
        """Splits the XML per variant and yields (id, XML) tuples."""
        soup = make_xml_soup(raw_response)
        for xml_element in soup.select(xml_element_tag):
            # FIXME: this will be extracted to a method:
            id_ = xml_element[xml_id_attribute]
            if db_name == 'snp':
                id_ = 'rs' + id_
            annotation = str(xml_element)
            yield (id_, annotation)

    def set_email_for_entrez(self):
        email_filepath = expanduser('~/.mail_address_for_Entrez')
        if not isfile(email_filepath):
            msg = ('Please set a mail for Entrez in {}. Entrez will notify '
                'that mail before banning you if your usage is too high.')
            raise Exception(msg)

        with open(email_filepath) as f:
            Entrez.email = f.read().strip()

