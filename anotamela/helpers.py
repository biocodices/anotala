import time
from os.path import isfile, expanduser
from itertools import zip_longest
import logging

from Bio import Entrez
from bs4 import BeautifulSoup
from tqdm import tqdm


logger = logging.getLogger(__name__)


def grouped(iterable, group_size):
    # Python recipe taken from:
    # https://docs.python.org/3.1/library/itertools.html#recipes
    args = [iter(iterable)] * group_size
    return ([e for e in t if e is not None] for t in zip_longest(*args))

def entrez_post_query(ids, db_name, batch_size, sleep_time, rettype='json',
                      xml_element_tag=None, xml_id_attribute=None):
    """
    Use Entrez POST query service to fetch a list of IDs in the given
    database. Yields tuples of (id, annotation).

    If rettype='xml', set xml_element_tag and xml_id_attribute to let the XML
    parser which XML elements are associated to each ID. For instance, for the
    db 'snp', set xml_element_tag='rs' and xml_id_attribute='rsid'. You have to
    manually explore the XML once in order to know those values.
    """
    if not Entrez.email:
        set_email_for_entrez()

    # Entrez POST queries are a two step process. You first POST your query
    # and get a WebEnv identifier and a QueryKey.
    logger.info('Create a serverside Job and get its ID')
    handle = Entrez.epost(db=db_name, id=','.join(ids))
    job = Entrez.read(handle)

    # Then you do a second query using those IDs in batches. More info:
    # http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec:entrez-webenv
    total = len(ids)
    msg = 'Fetch {} entries from Entrez DbSNP in batches of {}'
    logger.info(msg.format(total, batch_size))

    annotations = {}
    for offset in tqdm(list(range(0, total, batch_size))):
        # Fetch the raw response (XML or JSON)
        fetch_handle = Entrez.efetch(db=db_name, rettype=rettype,
                retstart=offset, retmax=batch_size,
                webenv=job['WebEnv'], query_key=job['QueryKey'])
        raw_response = fetch_handle.read()
        fetch_handle.close()

        if rettype == 'xml':
            # Split the XML per variant
            group_annotations = {}
            soup = make_xml_soup(raw_response)
            for xml_element in soup.select(xml_element_tag):
                # FIXME: this will be extracted to a method:
                id_ = xml_element[xml_id_attribute]
                if db_name == 'snp':
                    id_ = 'rs' + id_
                annotation = str(xml_element)
                yield (id_, annotation)

        #  elif rettype='json':
            #  for element in raw_response:
                ## extract ID!
                #  yield (id_, element)

        time.sleep(sleep_time)

def set_email_for_entrez():
    email_filepath = expanduser('~/.mail_address_for_Entrez')
    if not isfile(email_filepath):
        msg = ('Please set a mail for Entrez in {}. Entrez will notify '
               'that mail before banning you if your usage is too high.')
        raise Exception(msg)

    with open(email_filepath) as f:
        Entrez.email = f.read().strip()

def make_xml_soup(xml):
    return BeautifulSoup(xml, 'lxml')

