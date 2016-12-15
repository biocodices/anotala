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
    def _annotations_by_id(cls, _, pubmed_response):
        pubmed_records = pubmed_response['PubmedArticle'] + \
                         pubmed_response['PubmedBookArticle']
        yield from ((cls._extract_pmid(record), cls._parse_record(record))
                    for record in pubmed_records)

    @classmethod
    def _parse_annotation(cls, record):
        new_record = {}
        new_record['AMA_Citation'] = cls._generate_citation(record)
        new_record['CitationData'] = cls._generate_citation(record, as_dict=True)
        keys_to_keep = {
            'MedlineCitation.Article.Abstract.AbstractText': 'Abstract',
            'MedlineCitation.Article.ArticleDate': 'ArticleDate',
            'MedlineCitation.MeshHeadingList': 'Mesh',
            'PubmedData.ArticleIdList': 'Ids'
        }

        chosen_values = access_deep_keys(keys_to_keep.keys(), record,
                                         ignore_key_errors=True)
        for key, value in chosen_values.items():
            nicer_key = keys_to_keep[key]
            new_record[nicer_key] = value

        new_record['Ids'] = {i['IdType']: i['value']
                             for i in new_record['Ids']}

        if new_record['Mesh']:
            new_record['Mesh'] = [m['DescriptorName']['value']
                                  for m in new_record['Mesh']]

        return new_record

    @staticmethod
    def _extract_pmid(record):
        for id_ in record['PubmedData']['ArticleIdList']:
            if id_.attributes['IdType'] == 'pubmed':
                return str(id_)

    @staticmethod
    def _generate_citation(record, as_dict=False):
        article = record['MedlineCitation']['Article']

        citation_data = {
            'title': article['ArticleTitle'],
            'journal': article['Journal']['ISOAbbreviation'].replace('.', ''),
            'volume': article['Journal']['JournalIssue']['Volume'],
            'pages': article['Pagination']['MedlinePgn'],
        }

        try:
            citation_data['issue'] = article['Journal']['JournalIssue']['Issue']
        except KeyError:
            citation_data['issue'] = '-'

        try:
            citation_data['year'] = \
                article['Journal']['JournalIssue']['PubDate']['Year']
        except KeyError:
            citation_data['year'] = \
                article['Journal']['JournalIssue']['PubDate']['MedlineDate'].split()[0]

        if 'AuthorList' in article:
            author_list = []
            # Take the first three authors, not everybody
            for author in article['AuthorList'][:3]:
                if 'LastName' in author:
                    author_name = '{LastName} {Initials}'
                else:
                    author_name = '{CollectiveName}'
                author_list.append(author_name.format(**author))
            if len(article['AuthorList']) > 3:
                author_list.append('et al.')
            citation_data['authors'] = ', '.join(author_list)
        else:
            citation_data['authors'] = '[No authors listed]'

        tpl = '{authors}. {title}. {journal}. {year};{volume}({issue}):{pages}'
        citation = tpl.format(**citation_data).replace('..', '.')
        return citation_data if as_dict else citation

