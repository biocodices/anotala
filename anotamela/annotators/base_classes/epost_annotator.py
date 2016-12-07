from os.path import expanduser
import logging

from Bio import Entrez
from tqdm import tqdm

from anotamela.annotators.base_classes import AnnotatorWithCache


logger = logging.getLogger(__name__)
Entrez.tool = 'anotamela'


class EpostAnnotator(AnnotatorWithCache):
    """Base class for annotators that use Entrez epost service."""

    BATCH_SIZE = 1000

    def _batch_query_and_cache(self, ids):
        if hasattr(self, '_parse_id'):
            ids = [self._parse_id(id_) for id_ in ids]

        annotations = {}
        for batch_annotations in self._epost_batch_queries(ids):
            annotations.update(batch_annotations)
            self.cache.set(batch_annotations, namespace=self.SOURCE_NAME)

        return annotations

    def _epost_batch_queries(self, ids):
        """
        Use Entrez POST query service to fetch a list of IDs in the given
        database. Yields dictionaries with the annotations of each batch.

        If parse_xml=True, set xml_element_tag and xml_id_attribute to let the
        XML parser know which XML elements are associated to each ID. For
        instance, for the db 'snp', set xml_element_tag='rs' and
        xml_id_attribute='rsid'. You have to manually explore the XML once
        to know theese values.

        If self.USE_ENTREZ_READER is set, the raw response will be handled by
        Entrez.read() method, instead of returned as is.
        """
        if not Entrez.email:
            self.set_email_for_entrez()

        # Entrez POST queries are a two step process. You first POST your query
        # and get a WebEnv identifier and a QueryKey.
        logger.info('Create a serverside Job and get its ID')
        handle = Entrez.epost(db=self.ENTREZ_PARAMS['db'], id=','.join(ids))
        job_data = Entrez.read(handle)

        # Then you do a second query using those IDs, and you get the results
        # in batches. More info:
        # http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec:entrez-webenv
        total = len(ids)
        batch_size = min(self.BATCH_SIZE, total)
        msg = 'Fetch {} entries from Entrez DbSNP in batches of {}'
        logger.info(msg.format(total, min(batch_size, total)))
        for offset in tqdm(list(range(0, total, batch_size))):
            fetch_handle = Entrez.efetch(
                    db=self.ENTREZ_PARAMS['db'],
                    retmode=self.ENTREZ_PARAMS['retmode'],
                    webenv=job_data['WebEnv'], query_key=job_data['QueryKey'],
                    retstart=offset, retmax=batch_size
                )
            if hasattr(self, 'USE_ENTREZ_READER'):
                raw_response = Entrez.read(fetch_handle)
            else:
                raw_response = fetch_handle.read()
            fetch_handle.close()
            batch_annotations = self._annotations_by_id(raw_response)
            yield dict(batch_annotations)

    @staticmethod
    def set_email_for_entrez():
        email_filepath = expanduser('~/.mail_address_for_Entrez')

        try:
            with open(email_filepath) as f:
                Entrez.email = f.read().strip()
        except FileNotFoundError:
            msg = ('Please set a mail for Entrez in {}. Entrez will notify '
                   'that mail before banning you if your usage is too high.')
            raise FileNotFoundError(msg.format(email_filepath))

    @staticmethod
    def _annotations_by_id(raw_response):
        raise NotImplementedError

