import requests
import json
import time
import logging
import re

from tqdm import tqdm
import pandas as pd

from anotamela.annotators.base_classes import WebAnnotatorWithCache
from anotamela.helpers import grouped, path_to_source_file


logger = logging.getLogger(__name__)


class EnsemblAnnotator(WebAnnotatorWithCache):
    """
    Annotates rsids with Ensembl! REST service via POST requests.

    Set EnsemblAnnotator.full_info = True to get phenotypes, genotypes and
    population info besides the basic variant annotations.
    """
    SOURCE_NAME = 'ensembl'
    ANNOTATIONS_ARE_JSON = True

    BATCH_SIZE = 5
    # In theory Ensembl POST requests can handle up to 1,000 variants,
    # but every now and then I get timeouts and truncated responses when I try
    # to get the full info for as low as 25 variants at a time.
    # This is probably because we are asking for the full data: genotypes,
    # population genotypes, etc.

    SLEEP_TIME = 0

    api_version = 'GRCh37'
    full_info = False

    def _batch_query(self, ids):
        if self.proxies:
            logger.info('{} using proxies: {}'.format(self.name, self.proxies))

        for group_of_ids in tqdm(grouped(ids, self.BATCH_SIZE, as_list=True)):
            yield self._post_query(group_of_ids)
            time.sleep(self.SLEEP_TIME)

    def _post_query(self, ids):
        """
        Do a POST request to Ensembl REST api for a group of *ids*. Returns
        a dictionary with annotations per id. Requests should be done in
        batches of 1000 or less.
        """
        # No prefix needed for GRCh38
        url_prefix = 'grch37.' if self.api_version == 'GRCh37' else ''
        url = ('http://{}rest.ensembl.org/variation/homo_sapiens/?'
               .format(url_prefix))

        headers = {'Content-Type': 'application/json',
                   'Accept': 'application/json'}

        params = {'phenotypes': '1',
                  'genotypes': '1',
                  'pops': '1',
                  'population_genotypes': '1'}
        for key, value in params.items():
            url += '{}={};'.format(key, value)

        payload = {'ids': list(ids)}

        proxies = self.proxies or {}

        response = requests.post(url, headers=headers, proxies=proxies,
                                 data=json.dumps(payload))

        if response.ok:
            return response.json()
        else:
            logger.warn('Ensembl Error: {}'.format(response.text))
            response.raise_for_status()

    @classmethod
    def _parse_annotation(cls, annotation):
        if not cls.full_info:
            keys_to_remove = [
                'populations',
                'population_genotypes',
                'genotypes'
            ]
            for key in keys_to_remove:
                del(annotation[key])

        maf = annotation.get('MAF')
        if maf:
            annotation['MAF'] = float(maf)

        return annotation

    @classmethod
    def genotypes_list_to_dataframe(cls, genotypes, variant_id=None,
                                    ref_allele=None):
        """
        Given a list of Ensembl genotypes as annotated by EnsemblAnnotator
        with .full_info = True, return a tidy dataframe with one genotype per
        row, and extra columns for the sample name, population, superpopulation,
        a boolean indicating wether the sample belongs to 1000 Genomes or not,
        both alleles in a list.

        If marker_id is not None, the marker_id will also be set.

        If ref_allele is not None, the genotype will also be expressed as
        ALT allele dosage in a `alt_allele_dosage` column.

        1000 Genome compound sample identifiers like 1000GENOMES:phase_3:HG00097
        will be parsed.

        Example input:
            [{'gender': 'Female',
              'sample': '1000GENOMES:phase_3:HG00097',
              'genotype': 'C|T'},
               ...]

        """
        df = pd.DataFrame(genotypes)
        df['sample_full_name'] = df['sample']
        df['sample'] = df['sample'].map(cls.parse_1KG_sample_name)
        df['genotype_alleles'] = df['genotype'].map(cls.parse_genotype_string)
        df['in_1kg'] = df['sample_full_name'].str.contains('1000GENOMES')

        populations = cls.read_1kg_populations()
        df = df.merge(populations[['sample', 'pop', 'super_pop']],
                      how='left', on='sample')
        df = df.rename(columns={
            'pop': 'population',
            'super_pop': 'region'
        })

        df['alt_allele_dosage'] = None
        df['ref_allele'] = ref_allele
        if ref_allele:
            df['alt_allele_dosage'] = df['genotype_alleles'].map(
                lambda alleles: alleles.count(ref_allele)
            )
        df['variant_id'] = variant_id
        col_order = [
            'sample',
            'population',
            'region',
            'variant_id',
            'genotype_alleles',
            'ref_allele',
            'alt_allele_dosage',
            'in_1kg',
            'gender',
            'genotype',
            'sample_full_name',
        ]
        return df[col_order]

    @staticmethod
    def parse_1KG_sample_name(sample_name):
        if sample_name.startswith('1000GENOMES:'):
            sample_name = sample_name.split(':', 2)[2]
        return sample_name

    @staticmethod
    def parse_genotype_string(genotype_string):
        return tuple(re.split(r'\||/', genotype_string))

    @staticmethod
    def read_1kg_populations():
        fn = 'integrated_call_samples_v3.20130502.ALL.panel'
        df = pd.read_table(path_to_source_file(fn), sep='\s+')
        return df
