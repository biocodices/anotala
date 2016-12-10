import json
import logging

from anotamela.annotators.base_classes import EntrezAnnotator
from anotamela.helpers import access_deep_keys


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
    ENTREZ_PARAMS = {'db': 'pubmed', 'retmode': 'xml', 'service': 'epost'}
    USE_ENTREZ_READER = False

    @classmethod
    def _annotations_by_id(cls, _, pubmed_records):
        yield from ((cls._extract_pmid(record), cls._pythonify(record))
                    for record in pubmed_records)

    @classmethod
    def _parse_annotation(cls, record):
        new_record = {}
        new_record['AMA_Citation'] = cls._generate_citation(record)
        new_record['CitationData'] = cls._generate_citation(record, as_dict=True)
        keys_to_keep = {
            'MedlineCitation.Article.Abstract.AbstractText': 'Abstract',
            'MedlineCitation.Article.ArticleDate': 'ArticleDate',
            'MedlineCitation.Article.Journal': 'Journal',
            'PubmedData.ArticleIdList': 'Ids'
        }

        chosen_values = access_deep_keys(keys_to_keep.keys(), record,
                                         ignore_key_errors=True)
        for key, value in chosen_values.items():
            nicer_key = keys_to_keep[key]
            new_record[nicer_key] = value
        return new_record

    @staticmethod
    def _extract_pmid(record):
        for id_ in record['PubmedData']['ArticleIdList']:
            if id_.attributes['IdType'] == 'pubmed':
                return str(id_)

    @staticmethod
    def _pythonify(record):
        # FIXME: Awful hack. There must be a better way to do this.
        # The records are Bio.Entrez.Parser.StructureElement instances with
        # nested Bio.Entrez.Parser objects that would be a pain to manually
        # convert to their Python equivalents.
        # However, this serialization is non-optimal and it loses some data
        # saved in those elements' attributes. :(
        return json.loads(json.dumps(record))

    @staticmethod
    def _generate_citation(record, as_dict=False):
        article = record['MedlineCitation']['Article']
        author_list = ['{LastName} {Initials}'.format(**author)
                       for author in article['AuthorList'][:3]]
                       # ^ Take the first three authors, not everybody
        if len(article['AuthorList']) > 3:
            author_list.append('et al.')
        citation_data = {
            'authors': ', '.join(author_list),
            'title': article['ArticleTitle'],
            'journal': article['Journal']['ISOAbbreviation'].replace('.', ''),
            'year': article['Journal']['JournalIssue']['PubDate']['Year'],
            'volume': article['Journal']['JournalIssue']['Volume'],
            'issue': article['Journal']['JournalIssue']['Issue'],
            'pages': article['Pagination']['MedlinePgn'],
        }
        tpl = '{authors}. {title}. {journal}. {year};{volume}({issue}):{pages}'
        citation = tpl.format(**citation_data).replace('..', '.')
        return citation_data if as_dict else citation

