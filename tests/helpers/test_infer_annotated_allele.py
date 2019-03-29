import pytest

from anotala.helpers import infer_annotated_allele


@pytest.mark.parametrize('cds_change,expected_allele', [
    ('NM_1.1:c.123A>T', 'T'),
    ('NM_015001.2:c.*11199_*11199delTTT', 'delTTT'),
    ('chr20:g.62332640_62332641insG', 'insG'),
    ('c.123C=', 'C'),
    ('non_matching', 'non_matching'),
    ('', None),
])
def test_infer_annotated_allele(cds_change, expected_allele):
    assert infer_annotated_allele(cds_change) == expected_allele
