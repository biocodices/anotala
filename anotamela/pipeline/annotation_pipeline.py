import time
import logging
from itertools import chain
from os.path import expanduser
from functools import partial

import pandas as pd
import coloredlogs
from humanfriendly import format_timespan
from pprint import pformat

from anotamela.cache import create_cache, Cache
from anotamela.recipes import annotate_rsids_with_clinvar
from anotamela.annotators.base_classes.parallel_web_annotator import NoProxiesException
from anotamela.pipeline import (
    read_variants_from_vcf,
    annotate_rsids,
    extract_entrez_genes,
    annotate_entrez_gene_ids,
    get_omim_variants_from_entrez_genes,
    group_omim_variants_by_rsid,
    extract_pmids,
    annotate_pmids,
    update_pubmed_entries,
    extract_swissprot_ids,
    annotate_swissprot_ids,
    group_swissprot_variants_by_rsid,
    # annotate_clinvar_accessions,
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
        if isinstance(cache, Cache):
            cache = cache
        else:
            cache = create_cache(cache, **cache_kwargs)

        self.annotation_kwargs = {
            'cache': cache,
            'use_cache': use_cache,
            'use_web': use_web,
            'proxies': proxies,
            'sleep_time': sleep_time,
        }

        if proxies is None:
            raise NoProxiesException(
                "It's not advisable to run the complete pipeline "
                "without proxies, specially if you're going to "
                "annotate a lot of variants, because OMIM can get "
                "your IP banned. Try installing Tor locally and pass "
                "something this: "
                "proxies={'http': 'socks5://localhost:9050'}. If you "
                "still want to run without proxies, set proxies "
                "explicitely as an empty dict (proxies={})."
            )

    def run_from_vcf(self, vcf_path):
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
        self.vcf = vcf_path
        logger.info('Read "{}"'.format(self.vcf))

        variants = read_variants_from_vcf(expanduser(self.vcf))
        rs_variants = variants['rs_variants']
        other_variants = variants['other_variants']

        logger.info('{} variants with single rs'.format(len(rs_variants)))
        logger.info('{} other variants'.format(len(other_variants)))

        rs_annotations = self.run_from_rsids(rs_variants['id'].values)
        self.rs_variants = pd.merge(rs_variants, rs_annotations,
                                    left_on='id', right_on='id', how='left')

        self.other_variants = other_variants

    def run_from_rsids(self, rsids):
        """
        Given a list of rs IDs, annotate them and return the annotations in
        a pandas DataFrame. Variant annotations will be stored in
        self.rs_variants and related gene annotations in self.gene_annotations.
        """
        start_time = time.time()

        opts = pformat({**self.__dict__}, width=50)
        logger.info('Starting annotation pipeline with options:\n\n{}\n'
                    .format(opts))

        logger.info('Annotate the variants with rs ID')
        rs_variants = annotate_rsids(rsids, **self.annotation_kwargs)

        logger.info('Add ClinVar Variation Reports')
        clinvar_variations_per_rsid = annotate_rsids_with_clinvar(
            rsids,
            cache=self.annotation_kwargs['cache'],
            proxies=self.annotation_kwargs['proxies'],
            use_cache=self.annotation_kwargs['use_cache'],
            use_web=self.annotation_kwargs['use_web'],
            grouped_by_rsid=True,
        )
        rs_variants['clinvar_variations'] = \
            rs_variants['rsid'].map(clinvar_variations_per_rsid)

        logger.info('Extract Entrez gene data from the variants')
        dbsnp = rs_variants['dbsnp_myvariant']
        rs_variants['entrez_gene_ids'] = \
            dbsnp.fillna(False).apply(extract_entrez_genes, field='geneid')
        rs_variants['entrez_gene_symbols'] = \
            dbsnp.fillna(False).apply(extract_entrez_genes, field='symbol')

        logger.info('Annotate the Entrez genes associated to the variants')
        entrez_gene_ids = \
            list(chain.from_iterable(rs_variants['entrez_gene_ids']))
        gene_annotations = annotate_entrez_gene_ids(entrez_gene_ids,
                                                    **self.annotation_kwargs)

        logger.info('Get the OMIM variants described for the Entrez genes')
        omim_variants = \
            get_omim_variants_from_entrez_genes(entrez_gene_ids,
                                                **self.annotation_kwargs)

        logger.info('Extract PMIDs from the OMIM variants')
        pmids = extract_pmids(omim_variants)

        logger.info('Extract PMIDs from the GWAS Catalog entries')
        for gwas_entries in rs_variants['gwas_catalog'].dropna():
            pmids += extract_pmids(gwas_entries)

        logger.info('Annotate the PMIDs')
        pubmed_entries = annotate_pmids(pmids, **self.annotation_kwargs)

        logger.info("Update OMIM's PubMed entries with the PubMed annotations")
        omim_variants = update_pubmed_entries(omim_variants, pubmed_entries)

        logger.info('Associate each rs ID to a list of OMIM variants')
        rs_to_omim_variants = group_omim_variants_by_rsid(omim_variants)
        rs_variants['omim_entries'] = rs_variants['rsid'].map(rs_to_omim_variants)

        logger.info('Update GWAS PubMed entries with the PubMed annotations')
        func = partial(update_pubmed_entries,
                       pubmed_entries_by_pmid=pubmed_entries)
        rs_variants['gwas_catalog'] = \
            rs_variants['gwas_catalog'].map(func, na_action='ignore')

        logger.info('Extract Swissprot IDs')
        swissprot_ids = extract_swissprot_ids(gene_annotations['mygene'])

        logger.info('Get Swissprot variants from the swissprot gene IDs')
        swissprot_variants = annotate_swissprot_ids(swissprot_ids,
                                                    **self.annotation_kwargs)

        logger.info('Associate Swissprot variants to the rs IDs')
        rs_to_swissprot_variants = \
            group_swissprot_variants_by_rsid(swissprot_variants)
        rs_variants['uniprot_entries'] = \
            rs_variants['rsid'].map(rs_to_swissprot_variants)

        rs_variants.rename(columns={'clinvar_rs': 'clinvar_entries'},
                           inplace=True)

        self.rs_variants = rs_variants
        self.rs_variants.rename(columns={'rsid': 'id'}, inplace=True)

        self.gene_annotations = gene_annotations

        logger.info('Done! Took {} to complete the annotation'
                    .format(format_timespan(time.time() - start_time)))

        return self.rs_variants
