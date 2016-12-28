import pytest
from bs4 import BeautifulSoup

from anotamela import OmimGeneAnnotator


TEST_ENTRY_ID = '609708'
TEST_PHENO_ENTRY_ID = '100070'


@pytest.fixture(scope='module')
def annotator():
    annotator = OmimGeneAnnotator('mock_cache')
    annotator.PROXIES = {'http': 'socks5://beleriand.local:9150'}  # FIXME
    return annotator


@pytest.fixture(scope='module')
def html(annotator):
    # This will bring an OMIM page from the web, so that all the tests try
    # the current version of the HTML instead of a cached one
    annotation = annotator.annotate(TEST_ENTRY_ID, parse_data=False)
    return annotation[TEST_ENTRY_ID]


@pytest.fixture(scope='module')
def pheno_page_soup(annotator):
    html = annotator.annotate(TEST_PHENO_ENTRY_ID,
                              parse_data=False)[TEST_PHENO_ENTRY_ID]
    return BeautifulSoup(html, 'lxml')


@pytest.fixture(scope='module')
def soup(html):
    return BeautifulSoup(html, 'lxml')


@pytest.fixture(scope='module')
def variant_div(soup):
    # We'll test the first variant of the entry
    return soup.select('#allelicVariantsFold > div')[2:][0]  # FIXME: dupe code


def test_extract_entry_from_title(annotator, soup, pheno_page_soup):
    title_text = soup.title.text.strip()
    entry = annotator._extract_entry_from_title(title_text)
    assert entry['id'] == TEST_ENTRY_ID
    assert entry['type'] == 'gene_and_pheno'

    title_text = pheno_page_soup.title.text.strip()
    entry = annotator._extract_entry_from_title(title_text)
    assert entry['id'] == TEST_PHENO_ENTRY_ID
    assert entry['type'] == 'pheno_or_gene'


def test_extract_gene_from_soup(annotator, soup):
    gene = annotator._extract_gene_from_soup(soup)
    assert gene['symbol'] == 'LPL'
    assert gene['name'] == 'LIPOPROTEIN LIPASE'


def test_extract_variants_from_soup(annotator, soup, pheno_page_soup):
    variants = annotator._extract_variants_from_soup(soup)
    expected_total_variants = 42
    assert len(variants) == expected_total_variants

    expected_ids = ['{:04}'.format(n)
                    for n in range(1, expected_total_variants + 2)
                    if n != 22]  # Missing entry in our test ID
    assert [v['sub_id'] for v in variants] == expected_ids

    expected_keys = 'gene_omim_id gene_name gene_symbol'.split()
    assert all(variant[expected_key] for expected_key in expected_keys
                                     for variant in variants)

    # No variants in a Phenotype entry page:
    variants = annotator._extract_variants_from_soup(pheno_page_soup)
    assert not variants


def test_parse_variant_div(annotator, variant_div):
    entry = annotator._parse_variant_div(variant_div)

    expected_values = {
        'gene_symbol': 'LPL',
        'clinvar_accessions': ['RCV000001583'],
        'linked_mim_ids': ['238600'],
        'phenotype_names': ['LIPOPROTEIN LIPASE DEFICIENCY'],
        'prot_changes': ['ALA176THR'],
        'review_paragraphs': [
            'This mutation is sometimes called LPL-Bethesda.'  # Truncated
        ],
        'rsids': ['rs118204056'],
        'sub_id': '0001',
    }

    for key, expected_value in expected_values.items():
        if key == 'review_paragraphs':
            assert entry[key][0].startswith(expected_value[0])
        else:
            assert entry[key] == expected_value
