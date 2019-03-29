import re
from os.path import expanduser
from itertools import zip_longest

from Bio import Entrez

from anotala.helpers.variation_patterns import PROT_RE


def grouped(iterable, group_size, as_list=False):
    """Split an iterable in lists of <group_size> elements."""
    # Python recipe taken from:
    # https://docs.python.org/3.1/library/itertools.html#recipes
    args = [iter(iterable)] * group_size
    gen = ([e for e in t if e is not None] for t in zip_longest(*args))
    if as_list:
        return list(gen)
    else:
        return gen


def listify(maybe_list):
    if isinstance(maybe_list, list):
        return maybe_list
    return [maybe_list]


def parse_prot_change(prot_change):
    """
    Given a protein change like p.Ala123= or p.Ala123Ter, convert the '=' to
    the aminoacid ('Ala' in this case) or Ter to '*'.

    Leave the rest of the changes as they are.
    """
    match = PROT_RE.search(prot_change)

    if not match or re.search('(ins|del|dup)', prot_change):
        return prot_change

    info = match.groupdict()

    if info['aa2'] == '=':
        info['aa2'] = info['aa1']
    elif info['aa2'] == '*':
        info['aa2'] = 'Ter'
        # This unifies SnpEff's "*" with Clinvar, OMIM, and VEP use of "Ter"
        # for the stop codon.

    return 'p.{aa1}{pos}{aa2}'.format(**info)


def set_email_for_entrez():
    email_filepath = expanduser('~/.mail_address_for_Entrez')

    try:
        with open(email_filepath) as f:
            Entrez.email = f.read().strip()
    except FileNotFoundError:
        msg = ('Please set a mail for Entrez in {}. Entrez will notify '
               'that mail before banning you if your usage is too high.')
        raise FileNotFoundError(msg.format(email_filepath))


def camel_to_snake(s):
    """Convert a CamelCase string to a snake_case string."""
    return re.sub("([A-Z])", "_\\1", s).lower().lstrip('_')


def access_deep_keys(keys, dic, sep='.', ignore_key_errors=False):
    """
    Given a dictionary with maybe nested dictionaries in it, return a
    1-level flat new dictionary that includes the keys/values passed in keys.
    If a key has the sep character in it (e.g. 'Article.Abstract' if sep='.'),
    it will be treated as dic['Article']['Abstract'] and the new dicitonary
    will have a compound key new_dict['Article.Abstract'].

    Usage:

        > dic = {'foo': 1,
                 'bar': {'baz': 2,
                         'qux': {'ham': 3}}}
        > access_deep_keys(['foo', 'bar.baz', 'bar.qux.ham'], dic)
        # => {'bar.baz': 2, 'bar.qux.ham': 3, 'foo': 1}

    If ignore_key_errors=True, any failing keys will be ignored and you will
    get the result of the successful keys.
    """
    new_dic = {}
    for queried_key in keys:

        # A key with a separator character means we have to go deeper
        if sep in queried_key:
            key, deep_key = queried_key.split(sep, 1)
            if deep_key == '':
                msg = 'No key left when splitting "{}" with "{}"'
                raise ValueError(msg.format(queried_key, sep))

            if ignore_key_errors:
                deep_dic = dic.get(key, {})
            else:
                deep_dic = dic[key]

            if not isinstance(deep_dic, dict):
                msg = "You ask for '{}' but '{}' is not a dictionary"
                raise ValueError(msg.format(queried_key, key))

            deep_values = access_deep_keys([deep_key], deep_dic, sep=sep,
                                           ignore_key_errors=ignore_key_errors)
            new_dic[queried_key] = deep_values[deep_key]

        # A regular key with no separators in it
        else:
            if ignore_key_errors:
                new_dic[queried_key] = dic.get(queried_key)
            else:
                new_dic[queried_key] = dic[queried_key]

    return new_dic

