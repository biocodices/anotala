import pytest

from anotamela.helpers import (
    is_incidental_gene,
    is_incidental_pheno,
    listify,
)


INCIDENTAL_TEST_PARAMS = [
    # MIM ID, entry category, is incidental?
    (175100, 'pheno', True),
    (615373, 'pheno', False),
    (612048, 'gene', True),
    (605557, 'gene', False),
]


@pytest.mark.parametrize('mim_id,category,is_incidental', INCIDENTAL_TEST_PARAMS)
def test_is_incidental_gene(mim_id, category, is_incidental):
    check_func = {'pheno': is_incidental_pheno,
                  'gene': is_incidental_gene}[category]

    if is_incidental:
        assert check_func(mim_id)
        assert check_func(str(mim_id))
    else:
        assert not check_func(mim_id)
        assert not check_func(str(mim_id))


def test_listify():
    a_list = ['foo', 'bar']
    a_dict = {'foo': 'bar'}
    a_string = 'foo bar'

    assert listify(a_list) == a_list
    assert listify(a_dict) == [a_dict]
    assert listify(a_string) == [a_string]

