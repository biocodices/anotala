import re
import logging
from operator import itemgetter

from anotamela.annotators.base_classes import ParallelAnnotator
from anotamela.helpers import make_html_soup, gene_to_mim, mim_to_gene


logger = logging.getLogger(__name__)


class OmimGeneAnnotator(ParallelAnnotator):
    """
    Provider of OMIM annotations for genes (e.g. 605557).
    Responses are cached.

        > omim_gene_annotator = OmimGeneAnnotator()
        > omim_gene_annotator.annotate('605557)
        # => { '605557': ... }

    The base sleep time between requests is 60s because OMIM is quite strict
    against crawlers. However, you can lower this value at your own risk:

        > omim_variant_annotator.SLEEP_TIME = 30

    SLEEP_TIME is used as a base value that will be randomized each time.
    The user agent of each request is also randomized.
    """
    SOURCE_NAME = 'omim_genes'

    # OMIM is quite strict against web crawlers and it will ban your IP.
    # That's why we just avoid parallelization (BATCH_SIZE = 1) and we make
    # sure that there's a high random sleep time between requests:
    BATCH_SIZE = 1
    SLEEP_TIME = 60
    RANDOMIZE_SLEEP_TIME = True

    OMIM_URL = 'http://www.omim.org/entry/{}'
    PUBMED_URL = 'https://www.ncbi.nlm.nih.gov/pubmed/{}'
    REGEX = {
        'rsid': re.compile(r'dbSNP:(rs\d+)'),
        'prot_change': re.compile(r'\w+, (.+?)(?:,| \[| -| \()'),
        'aminoacids': re.compile(r'([A-Z]{3})-?(\d+)([A-Z]{3}|=)'),
    }

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

    def _query(self, mim_id):
        url = self._url(mim_id)
        return self._query_with_random_user_agent(url)

    @staticmethod
    def _url(mim_id):
        return 'http://omim.org/entry/{0}'.format(mim_id)

    @classmethod
    def _parse_annotation(cls, html):
        """Take the HTML from a OMIM Gene Entry and make a DataFrame with the
        variants for that gene."""

        variants = cls._extract_variants_from_html(html)
        phenotypes = cls._extract_phenotypes_from_html(html)
        references = cls._extract_references_from_html(html)

        for variant in variants:
            cls._add_data_to_variant(variant, phenotypes, references)

        return variants

    @classmethod
    def _add_data_to_variant(cls, variant, phenotypes, references):
        variant['variant_id'] = variant['gene_id'] + '.' + variant['sub_id']
        variant['rsid'] = '|'.join(cls.REGEX['rsid'].findall(variant['variant']))
        variant['entrez_id'] = mim_to_gene(variant['gene_id'])
        variant['gene_url'] = (cls.OMIM_URL.format(variant['gene_id']))

        matches = cls.REGEX['prot_change'].findall(variant['variant'])
        if matches:
            aminoacid_change = matches[0]

            aa_matches = cls.REGEX['aminoacids'].search(aminoacid_change)
            if aa_matches:
                # Change a protein change like GLY96ALA to Gly95Ala
                aa1, pos, aa2 = aa_matches.groups()
                aminoacid_change = 'p.{}{}{}'.format(aa1.capitalize(), pos,
                                                     aa2.capitalize())

            variant['prot_change'] = aminoacid_change

        variant['url'] = cls.OMIM_URL.format(
                variant['gene_id'] + '#' + variant['sub_id']
            )

        # Add the references found in the page to each variant,
        # if mentioned in the variant's review text.
        pmids = variant['pubmeds_summary'].values()
        pubmed_entries = [ref for ref in references
                          if ref.get('pmid') and ref.get('pmid') in pmids]
        for pubmed in pubmed_entries:
            pubmed['url'] = (cls.PUBMED_URL.format(pubmed['pmid']))
        variant['pubmeds'] = pubmed_entries

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

        for pheno in variant_phenotypes:
            pheno['url'] = cls.OMIM_URL.format(pheno['id'])

        # Hack to remove duplicated phenotypes
        unique_tuples = set(tuple(dic.items()) for dic in variant_phenotypes)
        unique_phenotypes = [dict(tup) for tup in unique_tuples]
        unique_phenotypes = sorted(unique_phenotypes, key=itemgetter('name'))

        variant['phenotypes'] = (unique_phenotypes or None)

    @classmethod
    def _extract_variants_from_html(cls, html):
        """
        Parses the HTML of an OMIM Gene Entry. Returns a list of variant
        entries found in the page.
        """
        html = html.replace('<br>', '<br/>')
        # ^ Need this so the parser doesn't think what comes after a <br>
        # is the <br>'s children. Added it to keep the newlines in OMIM
        # review texts.

        soup = make_html_soup(html)
        if soup.title.text.strip() == 'OMIM Error':
            return []

        omim_id = cls._extract_mim_id_from_soup(soup)
        if 'gene' not in omim_id['type']:
            return []

        gene = cls._extract_gene_from_soup(soup)
        variants = cls._extract_variants_from_soup(soup)

        for variant in variants:
            variant['gene_id'] = omim_id['id']
            variant['gene_name'] = gene['name']
            variant['gene_symbol'] = gene['symbol']

        return [variant for variant in variants if 'review' in variant]

    @classmethod
    def _extract_references_from_html(cls, html):
        soup = make_html_soup(html)
        references = []

        in_the_references_section = False
        for td in soup.select('td#floatingEntryContainer .wrapper-table td'):
            if not in_the_references_section:
                # Check again if the section started:
                in_the_references_section = (td.get('id') == 'references')
                continue

            # Get in the references section
            for td2 in td.select('.wrapper-table td.text'):
                if td2.get('id') and 'reference' in td2.get('id'):
                    continue

                fields = ['authors', 'title', 'publication', 'pubmed_str']
                values = [span.text.strip() for span in td2.select('span')]
                reference = dict(zip(fields, values))

                if 'pubmed_str' in reference:
                    pubmed_match = re.search(r'PubMed: (\d+)',
                                             reference['pubmed_str'])
                    if pubmed_match:
                        reference['pmid'] = pubmed_match.group(1)
                    del(reference['pubmed_str'])

                references.append(reference)

            break  # Leaving the references section, ignore the rest of <td>s

        return references

    @classmethod
    def _extract_phenotypes_from_html(cls, html):
        """
        Parse the HTML of an OMIM Entry page for a gene and return a list of
        phenotype dictionaries with name and id keys.
        """
        soup = make_html_soup(html)

        # Find the table title and advance to the next element
        table_title = soup.select('td#geneMap')

        if not table_title:  # Usually means we're in a Phenotype entry
            return []

        table_title = table_title[0]
        next_sib = [element for element in table_title.parent.next_siblings
                    if element != '\n'][0]

        phenotypes = []
        fields = 'name', 'id', 'inheritance', 'key'

        for tr in next_sib.select('.embedded-table tr')[1:]:  # Skips header
            row = [td.text.strip() for td in tr.select('td')
                   if 'rowspan' not in td.attrs]  # Skips Location column
            phenotype = dict(zip(fields, row))
            phenotypes.append(phenotype)

        return phenotypes

    @staticmethod
    def _extract_mim_id_from_soup(soup):
        mim_id = soup.select('span#title')[0].text.strip()
        # The MIM id is preceded by a * or #
        entry_symbol, entry_id = re.search(r'(\*|#|%|\+)?(\d+)', mim_id).groups()
        entry_types = {
            '*': 'gene',
            '#': 'phenotype',
            '%': 'pheno_or_gene',
            '+': 'gene_and_pheno',
            None: 'other',
        }
        return {'id': entry_id, 'type': entry_types[entry_symbol]}

    @staticmethod
    def _extract_gene_from_soup(soup):
        name_symbol = re.split(r';\s*', soup.select('td.title')[0].text.strip())
        if len(name_symbol) == 1:
            gene_name = name_symbol[0]
            # No symbol in the title --try in a different part of the page:
            hgcn_text = soup.select('td.subheading.italic.bold.text-font')
            gene_symbol = (hgcn_text and hgcn_text[0].text.split(': ')[-1])
        else:
            gene_name, gene_symbol = name_symbol

        return {'name': gene_name, 'symbol': gene_symbol}

    @staticmethod
    def _extract_variants_from_soup(soup):
        variants = []
        current_entry = {}

        in_the_variants_section = False
        table_rows = soup.select('td#floatingEntryContainer .wrapper-table td')
        for td in [td for td in table_rows if td.text]:
            inner_text = re.sub(r'\s+', ' ', td.text).strip()

            # Detect start and end of the ALLELIC VARIANTS section
            if not in_the_variants_section:
                in_the_variants_section = re.search(r'Table View', inner_text)
                continue
            elif re.search(r'REFERENCES', inner_text):
                variants.append(current_entry) if current_entry else None
                break

            # Title of new entry
            title_match = re.match(r'\.(\d{4}) (.+)', inner_text)
            if title_match:
                # Store the last entry
                variants.append(current_entry) if current_entry else None

                if re.search(r'REMOVED FROM|MOVED TO', inner_text):
                    current_entry = {}
                else:
                    current_entry = {
                            'sub_id': title_match.group(1),
                            'phenotype_names': [title_match.group(2)]
                    }
                continue

            if 'variant' not in current_entry:
                if '[ClinVar]' in inner_text:
                    current_entry['variant'] = inner_text.replace(' [ClinVar]',
                                                                  '')
                else:
                    phenos = [pheno for pheno in inner_text.split(', INCLUDED')
                              if pheno]
                    current_entry['phenotype_names'] += phenos
                continue

            # Rest of the entry will be the review

            # Extract the PubMed references from the text
            if 'pubmed' not in current_entry:
                current_entry['pubmeds_summary'] = {}

            for anchor in td.select('a.entry-reference'):
                current_entry['pubmeds_summary'][anchor.text] = \
                    anchor.get('pmid')

            # Extract the linked phenotype MIM IDs from the text
            if 'linked_mim_ids' not in current_entry:
                current_entry['linked_mim_ids'] = []

            omim_links = [anchor for anchor in td.select('a')
                          if anchor.get('href', '').startswith('/entry/')]
            for anchor in omim_links:
                current_entry['linked_mim_ids'].append(anchor.text)

            # Extract the linked phenotypes acronyms from the text
            if 'phenotype_symbols' not in current_entry:
                current_entry['phenotype_symbols'] = []

            pheno_symbols = re.findall(r'\((\w+); (?:see )?\d+\)', td.text)
            current_entry['phenotype_symbols'] += pheno_symbols

            # Get the review text itself without HTML elements
            def extract_review_texts(element):
                texts = []
                for child in element.children:
                    if child.name == 'div' and 'change' in child['class']:
                        # Nested div.change with the review texts
                        texts += extract_review_texts(child)
                    elif child.name == 'br':
                        texts.append('\n')
                    else:
                        try:
                            texts.append(child.text)
                        except AttributeError:
                            texts.append(child)

                return texts

            texts = ''.join(extract_review_texts(td))

            if 'review' not in current_entry:
                current_entry['review'] = texts
                continue

            current_entry['review'] += '\n\n' + texts

        return variants

