from itertools import zip_longest
from functools import lru_cache

from bs4 import BeautifulSoup


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

