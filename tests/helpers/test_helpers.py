import pytest
import types

from anotala.helpers import (
    is_incidental_gene,
    is_incidental_pheno,
    listify,
    parse_prot_change,
    grouped,
)


@pytest.mark.parametrize('prot_change,expected_result', [
    ('p.Ala123Gly', 'p.Ala123Gly'),
    ('p.Ala123=', 'p.Ala123Ala'),
    ('Ala123=', 'p.Ala123Ala'),
    ('p.Ala123*', 'p.Ala123Ter'),
    ('p.Ala123Ter', 'p.Ala123Ter'),
    ('p.Leu12_Leu14del', 'p.Leu12_Leu14del'),  # deletion is left untouched
    ('p.Ala160_Gln161insAla', 'p.Ala160_Gln161insAla'),  # insertion, same
    ('p.Pro56_Gly57insIleAlaPro', 'p.Pro56_Gly57insIleAlaPro'),
    ('p.Ala159dup', 'p.Ala159dup'),  # duplication, same
    ('Unknown', 'Unknown'),
])
def test_parse_prot_change(prot_change, expected_result):
    assert parse_prot_change(prot_change) == expected_result


@pytest.mark.parametrize('mim_id,category,is_incidental', [
    # MIM ID, entry category, is incidental?
    (175100, 'pheno', True),
    (615373, 'pheno', False),
    (612048, 'gene', True),
    (605557, 'gene', False),
])
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


def test_grouped():
    items = [1, 2, 3, 4]
    result = grouped(items, 2)
    expected_result = [[1, 2], [3, 4]]

    assert isinstance(result, types.GeneratorType)
    assert list(result) == expected_result

    result = grouped(items, 2, as_list=True)
    assert isinstance(result, list)
    assert result == expected_result

    result = grouped(items, 3, as_list=True)
    assert result == [[1, 2, 3], [4]]  # No None items, a shorter last group

