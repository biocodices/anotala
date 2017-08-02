from anotamela import ClinvarVariationAnnotator

import pytest


@pytest.fixture
def annotator():
    return ClinvarVariationAnnotator(cache='dict')


def test_annotate(annotator):
    result = annotator.annotate_one('1550')

    assert result['variation_id'] == '1550'
    assert 'c.953A>G' in result['variation_name']
    assert result['variation_type'] == 'Simple'

    assert len(result['genes']) == 1
    gene = result['genes'][0]
    assert gene['symbol'] == 'LPL'
    assert gene['entrez_id'] == '4023'
    assert gene['full_name'] == 'lipoprotein lipase'
    assert gene['strand'] == '+'
    assert gene['hgnc_id'] == '6677'

    assert len(result['clinical_assertions']) == 1
    clin = result['clinical_assertions'][0]
    assert clin['clinical_significance'] == 'Pathogenic'
    assert clin['date_last_submitted'] == '2010-12-30'
    assert clin['method'] == 'literature only'
    assert clin['submitter_name'] == 'OMIM'
    assert clin['phenotypes'] == [
        {'name': 'Hyperlipidemia, familial combined',
         'omim_id': '144250'}
    ]

