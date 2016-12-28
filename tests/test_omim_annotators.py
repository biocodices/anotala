from os.path import isfile

import pytest
from bs4 import BeautifulSoup

from anotamela import OmimGeneAnnotator
from helpers import get_test_file


TEST_ENTRY_ID = '609708'
TEST_PHENO_ENTRY_ID = '100070'


@pytest.fixture(scope='module')
def annotator(proxies):
    annotator = OmimGeneAnnotator('mock_cache')
    annotator.PROXIES = proxies
    return annotator


@pytest.fixture(scope='module')
def html(annotator):
    fn = get_test_file('htmls/omim_{}.html'.format(TEST_ENTRY_ID))
    with open(fn) as f:
        html = f.read()
    return html


@pytest.fixture(scope='module')
def pheno_page_soup(annotator):
    fn = get_test_file('htmls/omim_{}.html'.format(TEST_PHENO_ENTRY_ID))
    with open(fn) as f:
        html = f.read()
    return BeautifulSoup(html, 'lxml')


@pytest.fixture(scope='module')
def soup(html):
    return BeautifulSoup(html, 'lxml')


@pytest.fixture
def variant_divs(soup):
    return soup.select('#allelicVariantsFold > div')[2:]  # FIXME: dupe code


@pytest.fixture
def get_subdivs(annotator, variant_divs):
    # Hack to make a fixture that 'accepts arguments'

    def func(variant_index):
        return annotator._find_variant_subdivs(variant_divs[variant_index])

    return func


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


def test_parse_variant_div(annotator, variant_divs):
    # Test the first variant values
    entry = annotator._parse_variant_div(variant_divs[0])

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


def test_is_variant_ok(annotator, variant_divs):
    missing_variant_div = variant_divs[21]
    assert not annotator._is_variant_ok(missing_variant_div)

    present_variant_div = variant_divs[0]
    assert annotator._is_variant_ok(present_variant_div)


def test_find_variant_subdivs(annotator, get_subdivs):
    expected_sections = 'title extra description review empty'.split()
    assert all([section in get_subdivs(0) for section in expected_sections])


def test_parse_title_div(annotator, get_subdivs):
    title_div = get_subdivs(7)['title']
    result = annotator._parse_title_div(title_div)
    assert result['sub_id'] == '0008'
    assert result['phenotype_names'] == ['LIPOPROTEIN LIPASE DEFICIENCY']


def test_parse_extra_div(annotator, get_subdivs):
    extra_div = get_subdivs(18)['extra'][0]
    assert 'LPL-ARITA' in annotator._parse_extra_div(extra_div)


def test_parse_description_div(annotator, get_subdivs):
    desc_div = get_subdivs(0)['description']
    result = annotator._parse_description_div(desc_div)
    assert result['clinvar_accessions'] == ['RCV000001583']
    assert result['gene_symbol'] == 'LPL'
    assert result['prot_changes'] == ['ALA176THR']
    assert result['rsids'] == ['rs118204056']


def test_extract_review_paragraphs(annotator, get_subdivs):
    review_div = get_subdivs(0)['review']
    result = annotator._extract_review_paragraphs(review_div)
    assert result[0].startswith('This mutation is sometimes called')
    assert result[-1].endswith('as the chylomicronemia syndrome.')


def test_extract_pubmed_entries(annotator, get_subdivs):
    review_div = get_subdivs(1)['review']
    pubmeds = annotator._extract_pubmed_entries_from_review_div(review_div)
    expected_pmids = ['6645961', '1969408', '1975597', None, '2394828',
                      '1872917', '1351946', '11334614']
    assert [p['pmid'] for p in pubmeds] == expected_pmids


def test_extract_omim_links_from_review_div(annotator, get_subdivs):
    review_div = get_subdivs(0)['review']
    seen_mim_ids = annotator._extract_mim_ids_from_review_div(review_div)
    assert seen_mim_ids == ['238600']


def test_parse_empty_div(annotator, get_subdivs):
    empty_div = get_subdivs(0)['empty']
    annotator._parse_empty_div(empty_div)  # Expect no exceptions


def test_extract_phenotypes_from_soup(annotator, soup):
    phenos = annotator._extract_phenotypes_from_soup(soup)
    assert phenos == [
        {'id': '144250', 'inheritance': 'AD', 'key': '3',
         'name': 'Combined hyperlipidemia, familial'},
        {'id': '238600', 'inheritance': 'AR', 'key': '3',
         'name': 'Lipoprotein lipase deficiency'},
        {'id': '', 'inheritance': '', 'key': '3',
         'name': 'High density lipoprotein cholesterol level QTL 11'}
    ]


def test_extract_gene_from_soup(annotator, soup):
    gene = annotator._extract_gene_from_soup(soup)
    assert gene == {'name': 'LIPOPROTEIN LIPASE', 'symbol': 'LPL'}


@pytest.mark.skip(reason='TODO: write this test')
def test_add_data_to_variants(annotator):
    pass


def test_camelcase_prot_change(annotator):
    changes = {
        'ALA176THR': 'p.Ala176Thr',
        'INS': 'INS',
        'IVS2DS, G-A': 'IVS2DS, G-A',
        '6-KB DEL': '6-KB DEL',
        'GLN106TER': 'p.Gln106Ter',
    }
    for prot_change, expected_value in changes.items():
        assert annotator._camelcase_prot_change(prot_change) == expected_value
