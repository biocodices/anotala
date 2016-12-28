import pytest
from bs4 import BeautifulSoup

from anotamela import OmimGeneAnnotator


TEST_ENTRY = {
    'id': '605557',
    'type': 'gene',
    'gene_name': 'PR DOMAIN-CONTAINING PROTEIN 16',
    'gene_symbol': 'PRDM16',
}


@pytest.fixture(scope='session')
def annotator():
    annotator = OmimGeneAnnotator('mock_cache')
    annotator.PROXIES = {'http': 'socks5://beleriand.local:9150'}  # FIXME
    return annotator


@pytest.fixture(scope='session')
def html(annotator):
    # This will bring an OMIM page from the web, so that all the tests try
    # the current version of the HTML instead of a cached one
    annotation = annotator.annotate(TEST_ENTRY['id'], parse_data=False)
    return annotation[TEST_ENTRY['id']]


@pytest.fixture(scope='function')
def soup(html):
    return BeautifulSoup(html, 'html.parser')


def test_extract_entry_from_title(annotator, soup):
    title_text = soup.title.text.strip()
    entry = annotator._extract_entry_from_title(title_text)
    assert entry['id'] == TEST_ENTRY['id']
    assert entry['type'] == TEST_ENTRY['type']


def test_extract_gene_from_soup(annotator, soup):
    gene = annotator._extract_gene_from_soup(soup)
    assert gene['name'] == TEST_ENTRY['gene_name']
    assert gene['symbol'] == TEST_ENTRY['gene_symbol']


def test_extract_variants_from_soup(annotator, soup):
    variants = annotator._extract_variants_from_soup(soup)
    assert len(variants) == 6

    expected_ids = ['000{}'.format(n) for n in range(1, 7)]
    assert [v['sub_id'] for v in variants] == expected_ids

    expected_keys = 'gene_omim_id gene_name gene_symbol'.split()
    assert all(v[expected_key] for expected_key in expected_keys)

