import pytest
import pandas as pd

from anotamela.annotators import GeneEntrezLocalAnnotator

@pytest.fixture
def annotator():
    fp = pytest.helpers.file('Homo_sapiens.gene_info.mock.tsv')
    return GeneEntrezLocalAnnotator(
        path_to_annotations_file=fp
    )

def test_read_file(annotator):
    fp = pytest.helpers.file('Homo_sapiens.gene_info.mock.tsv')
    data = annotator._read_file(fp)

    assert data.loc[0]['GeneID'] == 1
    assert data.loc[0]['Symbol'] == 'GENE1'
    assert data.loc[1]['GeneID'] == 2
    assert data.loc[1]['Symbol'] == 'GENE2'


def test_parse_data(annotator):
    data = pd.DataFrame({'foo': [1, 2], 'bar': [2, 4]})
    result = annotator._parse_data(data)
    for ix, row in result.iterrows():
        original_row = data.loc[ix]
        for key in row.index:
            # Does no parsing!
            assert row[key] == original_row[key]


def test_annotate_many_ids(annotator):
    result = annotator._annotate_many_ids(['GENE1', '2', 1])
    assert len(result) == 2
    assert list(result['GeneID']) == [1, 2]
    assert list(result['Symbol']) == ['GENE1', 'GENE2']
