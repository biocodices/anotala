import re
import logging
from operator import itemgetter

from bs4 import BeautifulSoup

from anotamela.annotators.base_classes import ParallelWebAnnotator
from anotamela.helpers import (
    gene_to_mim,
    mim_to_gene,
    is_incidental_gene,
    is_incidental_pheno,
    parse_prot_change,
)


logger = logging.getLogger(__name__)


class OmimGeneAnnotator(ParallelWebAnnotator):
    """
    Provider of OMIM annotations for genes (e.g. 605557).
    Responses are cached.

        > omim_gene_annotator = OmimGeneAnnotator()
        > omim_gene_annotator.annotate_one('605557')
        # => { '605557': ... }

    The base sleep time between requests is 60s because OMIM is quite strict
    against crawlers. However, you can lower this value at your own risk:

        > omim_variant_annotator.SLEEP_TIME = 30

    SLEEP_TIME is used as a base value that will be randomized each time.
    The user agent of each request is also randomized.
    """
    SOURCE_NAME = 'omim_genes'

    BATCH_SIZE = 1
    # ^ OMIM guize are quite strict against web crawlers and they will ban your
    # IP. That's why we avoid parallelization (BATCH_SIZE=1) and we make
    # sure that there's a high random sleep time between requests:
    SLEEP_TIME = 20
    # However, OMIM *will ban your IP anyway* even with sleep times of around
    # 60 seconds, if you happen to leave this code scraping their site for a
    # couple of hours straight. So, be wise and don't use this without Tor
    # proxies!
    RANDOMIZE_SLEEP_TIME = True
    MANDATORY_PROXIES = True

    OMIM_URL = 'http://www.omim.org/entry/{}'
    PUBMED_URL = 'https://www.ncbi.nlm.nih.gov/pubmed/{}'
    PROT_RE = re.compile(r'([A-Z]{3})-?(\d+)([A-Z]{3}|=)')

    def annotate_from_entrez_ids(self, entrez_ids, **kwargs):
        entrez_ids = set(entrez_ids)
        omim_ids = [gene_to_mim(entrez_id) for entrez_id in entrez_ids
                    if entrez_id in gene_to_mim()]

        diff = len(entrez_ids) - len(omim_ids)
        if diff > 0:
            logger.warning('{} of the Entrez IDs have no MIM ID.'.format(diff))

        annotations = self.annotate(omim_ids, **kwargs)
        return {mim_to_gene(mim_id): annotation
                for mim_id, annotation in annotations.items()}

    @staticmethod
    def _url(mim_id):
        return 'http://omim.org/entry/{0}'.format(mim_id)

    @classmethod
    def _parse_annotation(cls, html):
        """Take the HTML from a OMIM Gene Entry and make a DataFrame with the
        variants for that gene."""

        # For some reason, 'lxml' works MUCH better than 'html.parser' for the
        # HTML found in OMIM pages:
        soup = BeautifulSoup(html, 'lxml')

        variants = cls._extract_variants_from_soup(soup)
        phenotypes = cls._extract_phenotypes_from_soup(soup)

        for variant in variants:
            cls._add_data_to_variant(variant, phenotypes)

        return variants

    @classmethod
    def _extract_variants_from_soup(cls, soup):
        """
        Extract all the allelic variants listed in the HTML soup of an OMIM
        *gene* entry page. Return them as dictionary records.
        """
        # Check it's a normal gene page and extract some general info:
        title_text = soup.title.text.strip()

        if title_text == 'OMIM Error':
            msg = 'Error parsing OMIM; page <title> says: "{}"'
            logger.warning(msg.format(title_text))
            return []

        omim_entry = cls._extract_entry_from_title(title_text)

        if 'gene' not in omim_entry['type']:
            msg = 'Not a gene OMIM entry; page <title> says: "{}"'
            logger.warning(msg.format(title_text))
            return []

        gene = cls._extract_gene_from_soup(soup)

        # Extract the variants from the page
        # First two <div>s of the alleles section are not variants, skip them
        variant_divs = soup.select('#allelicVariantsFold > div')[2:]
        variants = [cls._parse_variant_div(div) for div in variant_divs]
        variants = [v for v in variants if v and v['review_paragraphs']]

        for variant in variants:
            variant['gene_omim_id'] = omim_entry['id']
            variant['gene_name'] = gene['name']
            variant['gene_symbol'] = gene['symbol']

        return variants

    @staticmethod
    def _extract_entry_from_title(title_text):
        """
        Extract the MIM ID from the <title>, which is preceded by a symbol
        (*, #, %, +) that indicates what type of entry it is.

        For instance, '* 300644' is a gene entry with ID 30644.
        """
        title_re = re.compile(r'([%\*\+#])? (\d+)')
        # ^ Catches: '* 605557', '# 605557', '% 605557', '# 605557'
        entry_symbol, entry_id = title_re.search(title_text).groups()
        entry_types = {
            '*': 'gene',
            '#': 'phenotype',
            '%': 'pheno_or_gene',
            '+': 'gene_and_pheno',
            None: 'other',
        }
        return {'id': entry_id, 'type': entry_types[entry_symbol]}

    @classmethod
    def _parse_variant_div(cls, div):
        """
        Parse a <div> containing one allelic variant in an OMIM gene entry.
        """
        entry = {}

        if not cls._is_variant_ok(div):
            return entry

        subdivs = cls._find_variant_subdivs(div)

        title_info = cls._parse_title_div(subdivs['title'])
        entry.update(title_info)

        for extra_div in subdivs['extra']:
            if 'extra_names' not in entry:
                entry['extra_names'] = []
            entry['extra_names'] += cls._parse_extra_div(extra_div)

        description_info = cls._parse_description_div(subdivs['description'])
        entry.update(description_info)

        entry['review_paragraphs'] = \
            cls._extract_review_paragraphs(subdivs['review'])
        entry['pubmed_entries'] = \
            cls._extract_pubmed_entries_from_review_div(subdivs['review'])
        entry['linked_mim_ids'] = \
            cls._extract_mim_ids_from_review_div(subdivs['review'])

        cls._parse_empty_div(subdivs['empty'])

        return entry

    @staticmethod
    def _is_variant_ok(variant_div):
        title_div = variant_div.find('div')

        if re.search(r'(REMOVED FROM|MOVED TO)', title_div.text):
            return False

        # OMIM uses 'fake' ID 9999 for non-variant items
        if title_div.text.strip().startswith('.9999'):
            return False

        return True

    @staticmethod
    def _find_variant_subdivs(variant_div):
        """
        The typical variant has four <div>s:

            1- Title
            ?- [ optional extra <div>s with phenotype/variant names ]
            2- Shot description
            3- Review lengthy text
            4- One <div> that's empty

        Any number of extra <div>s might appear between title and description.
        Here we find the subdivs and identify them as the different sections
        based on their order. We return a dictionary like the following:

            {
                'title': <div> element,
                'review': <div> element,
                ...
             }
        """
        subdivs = variant_div.select('> div')
        categorized_subdivs = {}

        # The order of this popping is important, don't change it!
        categorized_subdivs['title'] = subdivs.pop(0)
        categorized_subdivs['empty'] = subdivs.pop()
        categorized_subdivs['review'] = subdivs.pop()
        categorized_subdivs['description'] = subdivs.pop()
        categorized_subdivs['extra'] = subdivs  # Remaining divs, possibly empty

        return categorized_subdivs

    @staticmethod
    def _parse_title_div(div):
        """Extract the allele ID (e.g. 0001, 0002, etc.) and usually the
        associated phenotype right next to it."""
        sub_id = div.select_one('span.mim-font').strong.text.strip()
        sub_id = sub_id[1:]  # Remove leading dot from IDs like ".0001"
        pheno_name = div.select_one('span.lookup').strong.text.strip()
        return {'sub_id': sub_id, 'phenotype_names': [pheno_name]}

    @staticmethod
    def _parse_extra_div(div):
        phenos = div.text.strip().split('\n')
        return [p.replace(', INCLUDED', '') for p in phenos]

    @staticmethod
    def _parse_description_div(div):
        """ Extract the gene symbol and aminoacid change from lines like
        'LPL, ASP204GLU'. Also get the dbSNP and RCV IDs from the links in the
        next line."""
        # Gene symbol and aminoacid change in the first line of text:
        description = div.text.strip().split('\n')[0].strip()
        gene_symbol, prot_changes = description.split(', ', maxsplit=1)
        prot_changes = prot_changes.split(' AND ')
        info = {'gene_symbol': gene_symbol, 'prot_changes': prot_changes}

        # dbSNP and ClinVar IDs + links
        for link in div.select('.btn-group > a'):
            if 'dbSNP' in link.text:
                if 'rsids' not in info:
                    info['rsids'] = []
                    info['ensembl_links'] = []
                info['rsids'].append(link.text.replace('dbSNP:', ''))
                info['ensembl_links'].append(link['href'])
            elif 'RCV' in link.text:
                assert 'clinvar_accessions' not in info
                # RCV accessions are all listed in the <a "title">
                # the link text only has the first one with ellipsis
                info['clinvar_accessions'] = link['title'].split(', ')
            elif 'ExAC' in link.text:
                if 'exac_ids' not in info:
                    info['exac_ids'] = []
                info['exac_ids'].append(link.text.replace('ExAC:', ''))
            else:
                msg = "I don't know how to parse this link: '{}'"
                raise NotImplementedError(msg.format(link.text))

        return info

    @staticmethod
    def _extract_review_paragraphs(div):
        return [p.text.strip() for p in div.select('p') if p.text.strip()]

    @staticmethod
    def _extract_pubmed_entries_from_review_div(div):
        pubmed_entries = []
        for link in div.select('a.entry-reference'):
            entry = {'pmid': link.get('pmid'),
                     'short_mention': link.text,
                     'citation': re.sub(r'\[.*?\]', '', link['title']).strip()}
            already_seen = any(pubmed['short_mention'] == entry['short_mention']
                               for pubmed in pubmed_entries)

            # Leave out entries with no PubMed ID
            if not already_seen and entry['pmid']:
                pubmed_entries.append(entry)

        return pubmed_entries

    @staticmethod
    def _extract_mim_ids_from_review_div(div):
        return [link.text for link in div.select('a')
                if link.get('href', '').startswith('/entry/')]

    @staticmethod
    def _parse_empty_div(div):
        assert not div.text.strip()

    @classmethod
    def _extract_phenotypes_from_soup(cls, soup):
        """Extract associated phenotypes info from the table with
        gene-phenotype relationships on top of the page."""
        phenotypes = []

        pheno_table = None
        for table in soup.select('table'):
            if table.text.strip().startswith('Location'):
                pheno_table = table
                break

        if not pheno_table:
            return phenotypes

        fields = 'name', 'id', 'inheritances', 'key'

        for tr in pheno_table.select('tbody tr'):  # Skips header
            row = [td.text.strip() for td in tr.select('td')
                   if 'rowspan' not in td.attrs]  # removes 'Location' column
            phenotype = dict(zip(fields, row))
            # Some phenotype names have surrounding curly braces or square
            # brackests (!) and might begin with a '?' sign (!!)
            phenotype['name'] = re.sub(r'[\[\]{}\?]', '', phenotype['name'])

            # Sometimes the inheritance comes with many comma-sep values
            # I have to tuple-ize the resulting list because otherwise
            # phenos can't be checked for uniqueness
            phenotype['inheritances'] = \
                tuple([inh for inh in phenotype['inheritances'].split(', ')
                       if inh])
            phenotypes.append(phenotype)

        return phenotypes

    @staticmethod
    def _extract_gene_from_soup(soup):
        """Extract the gene symbol from the <title> of the page or from
        some subheading in the body if that fails."""
        title_text = soup.title.text.strip()
        gene_regex = re.compile(r'\d+ - (.*)$')
        name_and_symbol = gene_regex.search(title_text).group(1).split('; ')
        if len(name_and_symbol) == 2:
            gene_name, gene_symbol = name_and_symbol
        elif len(name_and_symbol) == 1:
            # ^ No symbol in the title --try in a different part of the page:
            gene_name = name_and_symbol[0]
            hgcn_text = soup.select('td.subheading.italic.bold.text-font')
            gene_symbol = (hgcn_text and hgcn_text[0].text.split(': ')[-1])
        else:
            msg = "Couldn't parse this OMIM title for gene name + symbol: '{}'"
            raise Exception(msg.format(title_text))

        return {'name': gene_name, 'symbol': gene_symbol}

    @classmethod
    def _add_data_to_variant(cls, variant, phenotypes):
        variant['variant_id'] = variant['gene_omim_id'] + '.' + variant['sub_id']
        variant['gene_entrez_id'] = mim_to_gene(variant['gene_omim_id'])
        variant['gene_url'] = (cls.OMIM_URL.format(variant['gene_omim_id']))
        variant['incidental_gene'] = is_incidental_gene(variant['gene_omim_id'])
        variant['prot_changes'] = [cls._camelcase_prot_change(prot_change)
                                   for prot_change in variant['prot_changes']]
        variant['prot_changes'] = [parse_prot_change(prot_change)
                                   for prot_change in variant['prot_changes']]
        variant['url'] = cls.OMIM_URL.format(
                variant['gene_omim_id'] + '#' + variant['sub_id']
            )

        for pubmed in variant['pubmed_entries']:
            pubmed['url'] = (cls.PUBMED_URL.format(pubmed['pmid']))

        # Add the phenotypes cited in the table at the top of the page,
        # if they're mentioned in the variant's entry
        variant_phenotypes = [pheno for pheno in phenotypes
                              if pheno['id'] in variant['linked_mim_ids']]

        # This is a somewhat loose extra check of phenotype names
        # because the names are not always exactly the same
        # in the top table and in the variant entry.
        for pheno in phenotypes:
            for variant_pheno_name in variant['phenotype_names']:
                match1 = variant_pheno_name.lower() in pheno['name'].lower()
                match2 = pheno['name'].lower() in variant_pheno_name.lower()
                if match1 or match2:
                    variant_phenotypes.append(pheno)

        # Hack to remove duplicated phenotypes
        unique_tuples = set(tuple(dic.items()) for dic in variant_phenotypes)
        unique_phenotypes = [dict(tup) for tup in unique_tuples]
        unique_phenotypes = sorted(unique_phenotypes, key=itemgetter('name'))

        # Add extra info to the selected phenotypes
        # And revert a previous hack where the inheritance values where put
        # in a tuple instead of in a list
        for pheno in unique_phenotypes:
            pheno['url'] = cls.OMIM_URL.format(pheno['id'])
            pheno['incidental'] = is_incidental_pheno(pheno['id'])
            pheno['inheritances'] = list(pheno['inheritances'])  # Hack, see ^

        variant['phenotypes'] = (unique_phenotypes or [])

    @classmethod
    def _camelcase_prot_change(cls, prot_change):
        aa_matches = cls.PROT_RE.search(prot_change)
        if aa_matches:
            # Change a protein change like GLY96ALA to Gly96Ala
            aa1, pos, aa2 = aa_matches.groups()
            prot_change = 'p.{}{}{}'.format(aa1.capitalize(), pos,
                                            aa2.capitalize())
        return prot_change

