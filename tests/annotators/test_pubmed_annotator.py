from anotamela.annotators import PubmedAnnotator


def test_parse_annotation():
    record = {
        'MedlineCitation': {
            'Article': {
                'ArticleTitle': 'Title of the article',
                'ArticleDate': [{'Day': '01', 'Month': '01', 'Year': '1970'}],
                'AuthorList': [
                    {'CollectiveName': 'Some Consortium'},
                    {'LastName': 'Doe1', 'Initials': 'J'},
                    {'LastName': 'Doe2', 'Initials': 'J'},
                    {'LastName': 'Doe3', 'Initials': 'J'},
                    {'LastName': 'Doe4', 'Initials': 'J'},
                ],
                'Abstract': {'AbstractText': ['The text of the abstract.']},
                'Journal': {
                    'JournalIssue': {'Volumne': '1',
                                     'Issue': '1',
                                     'PubDate': {
                                         'Year': '1970',
                                         'Month': 'Jan',
                                     }},
                    'ISOAbbreviation': 'J. Clin. Invest.',
                },
                'Pagination': {'MedlinePgn': '100-20'},
            },
            'MeshHeadingList': [{'DescriptorName': {'value': 'mesh-term-1'}},
                                {'DescriptorName': {'value': 'mesh-term-2'}}],
        },
        'PubmedData': {
            'ArticleIdList': [
                {'IdType': 'pubmed', 'value': 'PMID-1'},
                {'IdType': 'pmc', 'value': 'PMC1'},
            ]
        }
    }

    result = PubmedAnnotator._parse_annotation(record)
    assert result['AMA_Citation'] == ('Some Consortium, Doe1 J, Doe2 J, et al. '
                                      'Title of the article. J Clin Invest. '
                                      '1970;(1):100-20')
    assert result['Abstract'] == ['The text of the abstract.']
    assert result['Mesh'] == ['mesh-term-1', 'mesh-term-2']
    assert result['pmid'] == 'PMID-1'

    assert 'url' in result
    assert 'ArticleDate' in result
    assert 'CitationData' in result

