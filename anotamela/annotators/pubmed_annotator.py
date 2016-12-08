import requests
from functools import lru_cache
import logging
from concurrent.futures import ThreadPoolExecutor

from anotamela.annotators.base_classes import EntrezAnnotator


logger = logging.getLogger(__name__)


class PubmedAnnotator(EntrezAnnotator):
    """
    Provider of Entrez PubMed summaries. Usage:

        > pubmed_entrez = PubmedAnnotator()
        > pubmed_entrez.annotate('23788249')  # Also accepts a list of IDs
        # => { '23788249': ... }
    """
    SOURCE_NAME = 'pubmed'
    ANNOTATIONS_ARE_JSON = True
    ENTREZ_PARAMS = {'db': 'pubmed', 'retmode': 'xml', 'service': 'esummary'}
    USE_ENTREZ_READER = True

    @staticmethod
    def _annotations_by_id(_, pubmed_records):
        yield from ((pubmed['Id'], pubmed) for pubmed in pubmed_records)

    @staticmethod
    def _parse_annotation(record):
        # Add AMA formatted citation
        authors = record['AuthorList'][:3]
        if len(record['AuthorList']) > 3:
            authors.append('et al.')
        record['AuthorsString'] = ', '.join(authors)
        template = '{AuthorsString}. {Title}. {Source}. {SO}'
        record['AMA_citation'] = template.format(**record).replace('..', '.')

        return record

