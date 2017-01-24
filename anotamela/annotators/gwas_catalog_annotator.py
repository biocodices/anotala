from itertools import groupby
from operator import itemgetter
import re

from anotamela.annotators.base_classes import ParallelAnnotator
from anotamela.helpers import camel_to_snake


class GwasCatalogAnnotator(ParallelAnnotator):
    SOURCE_NAME = 'gwas_catalog'
    ANNOTATIONS_ARE_JSON = True

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
            'betaNum',
            'shortForm',
            'catalogPublishDate',
            'betaDirection',
            'accessionId',
            'countriesOfRecruitment',
            'ancestralGroups',
            'ancestryLinks',
            'initialSampleDescription',
            'replicateSampleDescription',
            'reportedGene',
            'range',
            'title',
            'region',
            'chromLocation',
            'traitName_s',
            'traitUri',
            'parent',
            'chromosomeName',
            'riskFrequency',
        ]

    KEYS_TO_RENAME = {
            'rsids': 'rsId',
            'strongest_alleles': 'strongestAllele',
            'allele_impact': 'context',
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

    STRONGEST_ALLELE_REGEX = re.compile(r'(?P<rsid>rs\d+) ?-(?P<allele>.+)')

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

        associations = groups['association']
        return [cls._parse_association_entry(association)
                for association in associations]

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

        info = {k: v for k, v in info.items() if v}

        info['urls'] = {rsid: cls._web_url(rsid) for rsid in info['rsids']}

        if 'strongest_alleles' in info:
            info['genomic_alleles'] = []
            for rsid_allele in info['strongest_alleles']:
                info['genomic_alleles'] += cls._infer_allele(rsid_allele)

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

        info['pubmed_entries'] = cls._parse_pubmed_entries(info)

        return info

    @staticmethod
    def _web_url(rsid):
        """Return the GWAS Catalog web URL for an rsid."""
        return 'https://www.ebi.ac.uk/gwas/search?query={}'.format(rsid)

    @classmethod
    def _infer_allele(cls, rsid_alleles):
        """
        Given a 'strongest_allele' datum from GWAS Catalog like 'rs123-A' or
        'rs123-A; rs234-T', return a list of tuples with the rs IDs and the
        alleles: [('rs123', 'A'), ('rs234', 'T')].
        """
        rsid_alleles = rsid_alleles.split('; ')

        tuples = []
        for rsid_allele in rsid_alleles:
            match = cls.STRONGEST_ALLELE_REGEX.search(rsid_allele)
            if not match:
                raise ValueError("Couldn't parse the rsid_allele '{}'"
                                 .format(rsid_allele))
            allele = match.group('allele')
            tup = (match.group('rsid'), (None if allele == '?' else allele))
            tuples.append(tup)

        return tuples

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

    @staticmethod
    def _parse_ancestry_link(ancestry_link):
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
        na_values = ['NA', 'NR']
        for key in ['ancestries', 'countries', 'cities']:
            info[key] = [value for value in info[key]
                         if value and value not in na_values]

        # Include the raw input for checks of possibly missing data
        info['raw'] = ancestry_link

        return {k: v for k, v in info.items() if v and v not in na_values}

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

