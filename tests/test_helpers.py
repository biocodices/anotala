import pytest

from anotamela.helpers import is_incidental_gene, is_incidental_pheno


TEST_PARAMS = [
    (175100, 'pheno', True),
    (615373, 'pheno', False),
    (612048, 'gene', True),
    (605557, 'gene', False),
]

@pytest.mark.parametrize('mim_id,category,is_incidental', TEST_PARAMS)
def test_is_incidental_gene(mim_id, category, is_incidental):
    check_func = {'pheno': is_incidental_pheno,
                  'gene': is_incidental_gene}[category]

    if is_incidental:
        assert check_func(mim_id)
        assert check_func(str(mim_id))
    else:
        assert not check_func(mim_id)
        assert not check_func(str(mim_id))

