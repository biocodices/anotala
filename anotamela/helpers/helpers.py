import re
from os.path import expanduser
from itertools import zip_longest
from functools import lru_cache

from bs4 import BeautifulSoup
from Bio import Entrez


def grouped(iterable, group_size):
    # Python recipe taken from:
    # https://docs.python.org/3.1/library/itertools.html#recipes
    args = [iter(iterable)] * group_size
    return ([e for e in t if e is not None] for t in zip_longest(*args))


def make_xml_soup(xml):
    return BeautifulSoup(xml, 'lxml')


@lru_cache(maxsize=50)
def make_html_soup(html):
    return BeautifulSoup(html, 'html.parser')


def listify(maybe_list):
    if isinstance(maybe_list, list):
        return maybe_list
    return [maybe_list]


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

