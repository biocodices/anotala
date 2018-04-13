from itertools import groupby
from operator import itemgetter
import re

from anotamela.annotators.base_classes import ParallelWebAnnotator
from anotamela.helpers import camel_to_snake


class GwasCatalogAnnotator(ParallelWebAnnotator):
    SOURCE_NAME = 'gwas_catalog'
    ANNOTATIONS_ARE_JSON = True
    MANDATORY_PROXIES = True

    KEYS_TO_EXTRACT = [
            'id',
            'platform',
            'resourceName',
            'betaUnit',
            'resourcename',
            'reportedGeneLinks',
            'additionalAncestryDescription',
            'snpInteraction',
            'reportedGene',
            'publication',
            'publicationDate',
            'entrezMappedGeneLinks',
            'pValueExponent',
            'catalogPublishDate',
            'betaDirection',
            'accessionId',
            'countriesOfRecruitment',
            'ancestralGroups',
            'ancestryLinks',
            'initialSampleDescription',
            'replicateSampleDescription',
            'reportedGene',
            'title',
            'region',
            'traitUri',
            'parent',
        ]

    KEYS_TO_RENAME = {
            'rsids': 'rsId',
            'strongest_alleles': 'strongestAllele',
            'chrom_locations': 'chromLocation',
            'allele_impacts': 'context',
            'entrez_mapped_gene_symbols': 'entrezMappedGenes',
            'sample_size': 'numberOfIndividuals',
            'gwas_catalog_version': '_version_',
            'multisnp_haplotype': 'multiSnpHaplotype',
            'author_year_id': 'publicationLink',
            'pmid': 'pubmedId',
            'authors': 'author_s',
            'location_grch38': 'chromosomePosition',
            'p_value_base': 'pValueMantissa',
            'beta': 'betaNum',
            'ontology_ids': 'shortForm',
            'ontology_links': 'efoLink',  # needs parsing label|efo_id|efo_url
            'CI_range': 'range',
            'trait': 'traitName_s',
            'chromosome': 'chromosomeName',
            'risk_allele_frequency_in_controls': 'riskFrequency',
            'odds_ratio_per_copy': 'orPerCopyNum',
        }

    NA_VALUES = ['NA', 'NR']

    @staticmethod
    def _url(id_):
        # Gwas Catalog will take any search term here, but we typically
        # use it for rs IDs:
        url = 'https://www.ebi.ac.uk/gwas/api/search?q="{}"'.format(id_)

        # This grouping is optional, but facilitates the parsing afterwards:
        url += '&group=true&group.by=resourcename'
        return url

    @classmethod
    def _parse_annotation(cls, response):
        # This parsing depends on the URL above having the grouping options:
        groups = {group['groupValue']: group['doclist']['docs']
                  for group in response['grouped']['resourcename']['groups']}

        if not groups or 'association' not in groups:
            return

        return [cls._parse_association_entry(association)
                for association in groups['association']]

    @classmethod
    def _parse_association_entry(cls, association):
        """
        Given a single association entry from GWAS Catalog annotation,
        remove redundant information and parse the data. Returns a dictionary.
        """
        info = {}

        for key in cls.KEYS_TO_EXTRACT:
            new_key = camel_to_snake(key)
            info[new_key] = association.get(key)

        for new_key, key in cls.KEYS_TO_RENAME.items():
            info[new_key] = association.get(key)

        info = {k: v for k, v in info.items() if v and v not in cls.NA_VALUES}

        info['rsids'] = cls._parse_colon_separated_values(info['rsids'])

        # Some fields contain annotations for the different rs ids found in the
        # study. Split those and associate each annotation to the corresponding
        # rs ID:
        multiannotation_fields = [
            'strongest_alleles',
            'allele_impacts',
            'chrom_locations',
            'entrez_mapped_gene_symbols',
        ]
        for field in multiannotation_fields:
            if field not in info.keys():
                continue

            info[field] = cls._parse_colon_separated_values(info[field])
            info[field] = dict(zip(info['rsids'], info[field]))

        info['urls'] = {rsid: cls._web_url(rsid) for rsid in info['rsids']}

        float_key = 'risk_allele_frequency_in_controls'
        if float_key in info:
            info[float_key] = float(info[float_key])

        if 'strongest_alleles' in info:
            info['genomic_alleles'] = {}
            for rsid, rsid_allele in info['strongest_alleles'].items():
                info['genomic_alleles'][rsid] = cls._infer_allele(rsid_allele)

        if 'sample_size' in info:
            info['sample_size'] = cls._parse_sample_size(info['sample_size'])

        if 'entrez_mapped_gene_links' in info:
            mapped_genes = [cls._parse_entrez_mapped_gene_link(gene_link)
                            for gene_link in info['entrez_mapped_gene_links']]
            sorted_genes = \
                sorted(mapped_genes,
                       key=lambda gene: abs(gene['relative_position']))
            info['entrez_mapped_genes'] = sorted_genes
            del(info['entrez_mapped_gene_links'])

        if 'reported_gene_links' in info:
            reported_genes = [cls._parse_reported_gene_link(gene_link)
                              for gene_link in info['reported_gene_links']]
            info['reported_genes'] = reported_genes
            del(info['reported_gene_links'])

        if 'ancestry_links' in info:
            sample_info = [cls._parse_ancestry_link(ancestry_link)
                           for ancestry_link in info['ancestry_links']]
            info['sample_info'] = cls._group_sample_info(sample_info)
            del(info['ancestry_links'])

        if 'CI_range' in info:
            info['CI_range'] = cls._parse_ci_range(info['CI_range'])

        info['pubmed_entries'] = cls._parse_pubmed_entries(info)

        return info

    @staticmethod
    def _web_url(rsid):
        """Return the GWAS Catalog web URL for an rsid."""
        return 'https://www.ebi.ac.uk/gwas/search?query={}'.format(rsid)

    @classmethod
    def _infer_allele(cls, rsid_allele):
        """
        Given a 'strongest_allele' datum from GWAS Catalog like 'rs123-A',
        return a tuple with the rs IDs and the allele: ('rs123', 'A').
        """
        # Some IDs are given like rs1234NR, so we make that uniform with
        # the usual format, rs1234-NR:
        rsid_allele = rsid_allele.replace('NR', '-NR')
        pattern = r'(?P<rsid>rs\d+) ?-(?P<allele>.+)'
        match = re.compile(pattern).search(rsid_allele)

        if not match:
            print("Couldn't parse the rsid_allele '{}'".format(rsid_allele))
            return None

        allele = match.group('allele')

        return None if allele in ['?', 'NR'] else allele

    @staticmethod
    def _parse_pubmed_entries(association):
        """
        Given an association entry from GWAS Catalog with publication and pmid
        data, return a list of pubmed entries (each one a dictionary with a
        'pmid' key and more info on the entry if available).
        """
        pubmed_entries = {}
        for datum in association.get('author_year_id', []):
            author, year, pmid = [d.strip() for d in datum.split('|')]
            info = {'author': author,
                    'year': year,
                    'pmid': pmid}

            if pmid not in pubmed_entries:
                pubmed_entries[pmid] = info
            else:
                raise Exception('Repeated PubMed entry for this association {}'
                                .format(association))

        if 'pmid' in association:
            pmid = association['pmid']
            if pmid not in pubmed_entries:
                pubmed_entries[pmid] = {'pmid': pmid}

        return sorted(pubmed_entries.values(), key=itemgetter('pmid'))

    @staticmethod
    def _parse_sample_size(size_data):
        """
        Given the data of sample size for the association, return a dict with
        the same info explained. Example:

            > _parse_sample_size([100, 50])
              # => {'initial_study': 100, 'replication_study': 50}
        """
        return dict(zip(['initial_study', 'replication_study'], size_data))

    @staticmethod
    def _parse_entrez_mapped_gene_link(gene_link):
        """
        Given GWAS Catalog data of Entrez mapped gene link, return a
        dictionary with the same info explained. Example:

            > _parse_entrez_mapped_gene_links('HFE|3077|0|6')
              # => {'symbol': 'HFE',
                    'entrez_gene_id': '3077',
                    'relative_position': 0,
                    'chromosome': '6'}
        """
        fields = [
            'symbol',
            'entrez_gene_id',
            'relative_position',
            'chromosome',
        ]
        info = dict(zip(fields, gene_link.split('|')))
        info['relative_position'] = int(info['relative_position'])
        return info

    @staticmethod
    def _parse_reported_gene_link(gene_link):
        """
        Given GWAS Catalog data of Entrez reported gene link, return a
        dictionary with the same info explained. Example:

            > _parse_entrez_mapped_gene_links('HFE|3077|ENSG00000010704')
              # => {'symbol': 'HFE',
                    'entrez_gene_id': '3077',
                    'ensemble_id': 'ENSG00000010704'}
        """
        fields = [
            'symbol',
            'entrez_gene_id',
            'ensembl_id',
        ]
        return dict(zip(fields, gene_link.split('|')))

    @classmethod
    def _parse_ancestry_link(cls, ancestry_link):
        """
        Given an ancestry link from GWAS Catalog, return a dictionary with the
        parsed information. Example:

            > link = ('initial|Netherlands|Netherlands|European|911|Utrecht, '
                      'Netherlands; St Radboud, Netherlands;')
            > _parse_ancestry_link(link)
            > # => {'study': 'initial',
                    'extra': 'Netherlands',
                    'countries': ['Netherlands'],
                    'ancestries': ['European'],
                    'sample_size': 911,
                    'cities': ['Utrecht, Netherlands',
                               'St RAdboud, Netherlands']}
        """
        fields = [
            'study_type',
            'extra',
            'countries',
            'ancestries',
            'sample_size',
            'cities',
        ]
        info = dict(zip(fields, ancestry_link.split('|')))
        info['sample_size'] = int(info['sample_size'])
        info['ancestries'] = info['ancestries'].split(',')
        info['countries'] = info['countries'].split(',')
        info['cities'] = [city.strip() for city in info['cities'].split(';')]

        # Remove NA values inside lists
        for key in ['ancestries', 'countries', 'cities']:
            info[key] = [value for value in info[key]
                         if value and value not in cls.NA_VALUES]

        # Include the raw input for checks of possibly missing data
        info['raw'] = ancestry_link

        return {k: v for k, v in info.items() if v and v not in cls.NA_VALUES}

    @staticmethod
    def _group_sample_info(sample_info):
        """
        Given samples information parsed by _parse_ancestry_link, groups the
        entries by study type (i.e. 'initial' vs 'replication'). Returns a
        dictionary.
        """
        return {study_type: list(studies)
                for study_type, studies
                in groupby(sample_info, itemgetter('study_type'))}

    @staticmethod
    def _parse_colon_separated_values(values):
        """
        Given a list of values separated by colons, as given by GWAS Catalog
        annotation, parse them and return the values in the same order.
        Returns a list.

        Example:

            > _parse_colon_separated_values(['intron_variant; intron_variant'])
              # => ['intron_variant', 'intron_variant']
            > _parse_colon_separated_values(['rs123; rs234'])
              # => ['rs123', 'rs234']
        """
        parsed_values = []

        for value in values:
            # some are simple like 'intron_variant', others are compound like
            # 'intron_variant; intron_variant; 3_prime_UTR_variant',
            # so we need to split those:
            parsed_values += value.split('; ')

        return parsed_values

    @classmethod
    def _parse_ci_range(cls, ci_range):
        """
        Given a confidence interval range as a string like '[1.01-1.25]',
        parse it to get the numbers in a dict like:
        {'lower_limit': 1.01, 'upper_limit': 1.25}.

        Return None for '[NR]' values.
        """
        if ci_range == '[NR]':
            return None

        # Formats for this range are inconsistent. So far I've seen:
        #
        # [1.01-1.25]
        # 1.01-1.25
        # 1.01 - 1.25
        # [1.72-4]
        #
        # So I just try to capture the two numbers:
        parsed = ci_range.replace('[', '').replace(']', '')
        numbers = [n.strip() for n in parsed.split('-')]
        numbers = [float(n) for n in numbers]
        values = dict(zip(['lower_limit', 'upper_limit'], numbers))

        if not len(values) == 2:
            msg = "Clouldn't find two confidence interval values in \"{}\""
            raise ValueError(msg.format(ci_range))

        return values

