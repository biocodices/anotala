import pytest

import pandas as pd
import numpy as np

from anotamela.pipeline import (
    fix_genomic_alleles_for_variant,
    fix_genomic_allele_given_VCF_alleles
)


def test_fix_genomic_alleles_for_variant():
    variant = pd.Series({
        'ref': 'A',
        'alt': ['AC'],
        'string': 'string',
        'int': 1,
        'entries': [{'genomic_allele': 'insC'}],
        'entry': {'genomic_allele': 'insC'},
    })
    result = fix_genomic_alleles_for_variant(variant)
    assert result['string'] == 'string'
    assert result['int'] == 1
    assert result['entries'] == [{'genomic_allele': 'AC'}]
    assert result['entry'] == {'genomic_allele': 'AC'}


def test_fix_genomic_allele_given_VCF_alleles():
    f = fix_genomic_allele_given_VCF_alleles

    insertion = {'genomic_allele': 'insC'}
    result = f(insertion, ref='A', alts=['AC'])
    assert result['genomic_allele'] == 'AC'

    insertion = {'genomic_allele': 'insCC'}
    result = f(insertion, ref='A', alts=['AC', 'ACC', 'ACCC'])
    assert result['genomic_allele'] == 'ACC'

    deletion = {'genomic_allele': 'del'}
    result = f(deletion, ref='AC', alts=['A'])
    assert result['genomic_allele'] == 'A'

    result = f(deletion, ref='ACC', alts=['A', 'AC'])
    assert result['genomic_allele'] == 'del' # If ambiguous, don't choose

    multiple_entries = [{'genomic_allele': 'insC'},
                        {'genomic_allele': 'insCC'}]
    result = f(multiple_entries, ref='A', alts=['AC', 'ACC'])
    assert result[0]['genomic_allele'] == 'AC'
    assert result[1]['genomic_allele'] == 'ACC'

    with pytest.raises(ValueError):
        f('a string', ref='A', alts=['AC'])

    with pytest.raises(ValueError):
        f(['a', 'list', 'of', 'strings'], ref='A', alts=['AC'])

    # Doesn't fix missing genomic allele, and leaves other values untouched:
    no_genomic_allele = {'other_key': 'other_value'}
    result = f(no_genomic_allele, ref='A', alts=['AC'])
    assert result == no_genomic_allele

    # Deal with NaN
    result = f(np.nan, ref='A', alts=['AC'])
    assert result is None
