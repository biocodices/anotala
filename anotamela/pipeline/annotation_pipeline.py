import time
import logging
from itertools import chain
from functools import partial
from os.path import expanduser

import pandas as pd
import coloredlogs
from humanfriendly import format_timespan
from pprint import pformat

from anotamela.cache import create_cache, Cache
from anotamela.annotators import RS_ANNOTATOR_CLASSES
from anotamela.pipeline import (
    read_variants_from_vcf,
    annotate_rsids,
    extract_entrez_genes,
    annotate_entrez_gene_ids,
)


logger = logging.getLogger(__name__)
coloredlogs.DEFAULT_LOG_FORMAT = '[@%(hostname)s %(asctime)s] %(message)s'
coloredlogs.install(level='INFO')


class AnnotationPipeline:
    def __init__(self, cache, use_cache=True, use_web=True, proxies=None,
                 sleep_time=None, **cache_kwargs):
        """
        Initialize a pipeline with a given set of options. The options will
        be used for any subsequent pipeline.run() actions.

        - cache is mandatory. You can 'manually' instantiate a Cache
          (either PostgresCache or RedisCache) and pass it here, or you can
          specify 'redis' or 'postgres' and let the pipeline do that for you.
        - use_cache (default=True): whether to use or not data found in cache
          for each variant.
        - use_web (default=True): whether to use or not web data to annotate
          the variants. If use_cache is also set, the web will be used only
          to annotate the ones not found in cache. If use_cache=False, every
          variant will be annotated from web, updating any previous cached
          data for those variants.
        - proxies (default=None) is optional. If set, it should be a dictionary
          of proxies that will be used by the requests library. For instance:
          {'http': 'socks5://localhost:9050'}
        - sleep_time (default=None) is optional. If set, it will be used to
          override all annotators SLEEP_TIME between queries.
        - **cache_kwargs will be passed to the Cache constructor if the cache
          option is not already a Cache instance.

        See the docstring of Pipeline.run for some usage examples.
        """
        self.use_cache = use_cache
        self.use_web = use_web
        self.proxies = proxies
        self.sleep_time = sleep_time

        if isinstance(cache, Cache):
            self.cache = cache
        else:
            self.cache = create_cache(cache, **cache_kwargs)

        self.annotation_kwargs = {
            'cache': self.cache,
            'use_cache': self.use_cache,
            'use_web': self.use_web,
            'proxies': self.proxies,
            'sleep_time': self.sleep_time,
        }

    def run(self, vcf_path):
        """
        Annotate the given VCF file (accepts gzipped files too). Returns
        a DataFrame of annotations per rs ID. The IDs that weren't regular
        rs IDs and thus were not annotated are stored in self.other_variants

        Examples:

            > from anotamela import Pipeline

            > # Annotate using only cached variants
            > pipeline = Pipeline(cache='postgres', use_web=False)
            > pipeline.run('~/variants.vcf.gz')

            > # Or use web too, and pass some args for the cache:
            > pipeline = Pipeline(cache='postgres',
                                  credentials_filepath='~/.pg_creds.yml')
            > pipeline.run('~/variants.vcf.gz')

            > # Use Redis for cache:
            > pipeline = Pipeline(cache='redis', host='192.168.1.10')
            > # ^ Don't pass the host argument to use 'localhost'
            > pipeline.run('~/variants.vcf')

            # use_cache=False will reannotate any already cached variants:
            > pipeline = Pipeline(cache='postgres', use_cache=False)
            > pipeline.run('~/variants.vcf')

        """
        self.start_time = time.time()

        opts = pformat({**self.__dict__, 'vcf_path': vcf_path} , width=50)
        logger.info('Starting annotation pipeline with options:\n\n{}\n'
                    .format(opts))

        logger.info('Read "{}"'.format(vcf_path))

        variants = read_variants_from_vcf(expanduser(vcf_path))
        self.rs_variants = variants['rs_variants']
        self.other_variants = variants['other_variants']

        logger.info('{} variants with single rs'.format(len(self.rs_variants)))
        logger.info('{} other variants'.format(len(self.other_variants)))

        logger.info('Annotate the variants with rs ID')
        rs_annotations = annotate_rsids(self.rs_variants['id'],
                                        **self.annotation_kwargs)
        self.rs_variants = pd.merge(self.rs_variants, rs_annotations,
                                    on='id', how='left')

        logger.info('Extract Entrez gene data from the variants')

        dbsnp = self.rs_variants['dbsnp_myvariant']
        self.rs_variants['entrez_gene_ids'] = \
                dbsnp.fillna(False).apply(extract_entrez_gene, field='geneid')
        self.rs_variants['entez_gene_symbols'] = \
                dbsnp.fillna(False).apply(extract_entrez_gene, field='symbol')

        logger.info('Annotate the Entrez genes associated to the variants')

        entrez_gene_ids = chain.from_iterable(self.rs_variants['entrez_gene_ids'])
        self.gene_annotations = annotate_entrez_genes(list(entrez_gene_ids),
                                                      **self.annotation_kwargs)

        self._add_omim_variants_info()
        self._add_pubmed_entries_to_omim_entries()
        self._add_uniprot_variants_info()

        self.rs_variants.rename(columns={'clinvar_rs': 'clinvar_entries'},
                                inplace=True)

        self.end_time = time.time()
        elapsed = format_timespan(self.end_time - self.start_time)
        logger.info('Done! Took {} to complete the pipeline'.format(elapsed))

        return self.rs_variants

    def _add_omim_variants_info(self):
        # Get all OMIM variants in the affected genes
        annotator = OmimGeneAnnotator(cache=self.cache)
        annotator.PROXIES = self.proxies
        omim_variants_per_gene = annotator.annotate_from_entrez_ids(
                self.gene_annotations['entrez_id'],
                use_cache=self.use_cache,
                use_web=self.use_web,
            )

        # Keep the variants with the rs IDs that we are interested in
        rs_to_omim_variants = {}

        omim_variants = [v for gene_variants in omim_variants_per_gene.values()
                           for v in gene_variants]

        for omim_variant in omim_variants:
            # One OMIM variant might correspond to many
            # rs IDs, so we need to check each rs associated to an OMIM
            # variant.
            for rs in omim_variant.get('rsids', []):
                if rs not in self.rs_variants['id'].values:
                    continue

                if rs not in rs_to_omim_variants:
                    # There might be more than one OMIM variant for a given rs
                    # --the typical case is the one of multiallelic SNPs.
                    # Hence, we return the OMIM annotation for an rs always as
                    # a *list* of entries.
                    rs_to_omim_variants[rs] = []

                rs_to_omim_variants[rs].append(omim_variant)

        self.rs_variants['omim_entries'] = \
                self.rs_variants['id'].map(rs_to_omim_variants)

    def _add_pubmed_entries_to_omim_entries(self):
        pubmed_ids = set()

        for omim_entries in self.rs_variants['omim_entries'].dropna():
            for entry in omim_entries:
                for pubmed in entry.get('pubmed_entries', []):
                    pubmed_ids.add(pubmed['pmid'])

        pubmed_annotations = \
                PubmedAnnotator(cache=self.cache).annotate(pubmed_ids)

        for omim_entries in self.rs_variants['omim_entries'].dropna():
            for entry in omim_entries:
                for pubmed in entry.get('pubmed_entries', []):
                    extra_info = pubmed_annotations[pubmed['pmid']]
                    pubmed.update(extra_info)

    def _add_uniprot_variants_info(self):
        uniprot_ids = set()

        for gene in self.gene_annotations['mygene'].fillna(False):
            if gene and 'swissprot' in gene:
                ids = gene['swissprot']
                if not isinstance(ids, list):
                    ids = [ids]
                uniprot_ids = uniprot_ids | set(ids)

        annotator = UniprotAnnotator(cache=self.cache)
        annotator.PROXIES = self.proxies
        uniprot_annotations = annotator.annotate(
            uniprot_ids,
            use_cache=self.use_cache,
            use_web=self.use_web
        )
        uniprot_variants = list(chain.from_iterable(uniprot_annotations.values()))

        if uniprot_variants:
            uniprot_df = pd.DataFrame(uniprot_variants)
            uniprot_df = uniprot_df.set_index('rsid', drop=False)

            def rs_to_uniprot_variants(rs):
                if rs not in uniprot_df.index:
                    return None
                return uniprot_df.loc[[rs]].to_dict('records')

            self.rs_variants['uniprot_entries'] = \
                self.rs_variants['id'].map(rs_to_uniprot_variants)

        else:
            logger.warning('No uniprot variants annotated?')

