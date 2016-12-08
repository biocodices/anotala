import logging

from Bio import Entrez
from tqdm import tqdm

from anotamela.annotators.base_classes import AnnotatorWithCache
from anotamela.helpers import set_email_for_entrez


logger = logging.getLogger(__name__)


class EntrezAnnotator(AnnotatorWithCache):
    """
    Base class for annotators that use one of Entrez services. Classes that
    inherit from this one should have:

        - class variable SOURCE_NAME

        - class variable ENTREZ_PARAMS =
            { 'db': ..., 'retmode': 'json' , 'service': (epost|esummary) }

            'service' lets you choose between epost and esummary services.

        - a method _annotations_by_id(raw_response) that takes the raw
          response with info for many ids and yields tuples of (id, data) for
          each of the queried IDs, thus splitting the response per ID.

        - optionally, a _parse_id(id_) method that will transform the IDs in
          any way prior to the query. This is useful, for instance, to remove
          'rs' from rs IDs.

        - optionally, a USE_ENTREZ_READER class variable to indicate that the
          response from Entrez should be handled by Entrez.read(). This works
          for some but not all DBs.
    """
    def _batch_query(self, ids):
        """
        Use Entrez query service to fetch a list of IDs in batches, in the
        database defined in self.ENTREZ_PARAMS['db']. Yields dictionaries with
        the annotations of each batch.

        If self.USE_ENTREZ_READER is set, the raw response will be handled by
        Entrez.read() method, instead of returned as is.
        """
        if hasattr(self, '_parse_id'):
            ids = [self._parse_id(id_) for id_ in ids]

        if not Entrez.email:
            set_email_for_entrez()

        total = len(ids)
        self.batch_size = 1000 if self.ENTREZ_PARAMS['service'] == 'epost' else 200
        msg = 'Fetch {} entries from Entrez "{}" in batches of {}'.format(
            total, self.ENTREZ_PARAMS['db'], min(self.batch_size, total))
        logger.info(msg)

        if self.ENTREZ_PARAMS['service'] == 'epost':
            handles = self._epost_query(ids)
        elif self.ENTREZ_PARAMS['service'] == 'esummary':
            handles = self._esummary_query(ids)
        else:
            msg = "I don't have a query method for '{}'"
            raise NotImplementedError(msg.format(self.ENTREZ_PARAMS['service']))

        for handle in handles:
            if hasattr(self, 'USE_ENTREZ_READER'):
                raw_response = Entrez.read(handle)
            else:
                raw_response = handle.read()

            handle.close()
            yield from self._annotations_by_id(ids, raw_response)

    def _esummary_query(self, ids):
        handle = Entrez.esummary(db=self.ENTREZ_PARAMS['db'], id=','.join(ids))
        yield handle

    def _epost_query(self, ids):
        # Entrez POST queries are a two step process. You first POST your query
        # and get a WebEnv identifier and a QueryKey.
        logger.info('Create a serverside Job and get its ID')
        handle = Entrez.epost(db=self.ENTREZ_PARAMS['db'], id=','.join(ids))
        job_data = Entrez.read(handle)

        # Then you do a second query using the job data, and you get the
        # results in batches. More info:
        # http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec:entrez-webenv
        logger.info('Get the results from the job')
        for offset in tqdm(list(range(0, len(ids), self.batch_size))):
            fetch_handle = Entrez.efetch(
                    db=self.ENTREZ_PARAMS['db'],
                    retmode=self.ENTREZ_PARAMS['retmode'],
                    webenv=job_data['WebEnv'], query_key=job_data['QueryKey'],
                    retstart=offset, retmax=self.batch_size
                )

            yield fetch_handle

    @staticmethod
    def _annotations_by_id(ids, raw_response):
        raise NotImplementedError

