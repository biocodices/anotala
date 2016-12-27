import time
import logging
from itertools import chain

import pandas as pd
import coloredlogs
from vcf_to_dataframe import vcf_to_dataframe
from humanfriendly import format_timespan
from pprint import pformat

from anotamela.cache import create_cache, Cache
from anotamela.annotators import *


logger = logging.getLogger(__name__)
coloredlogs.DEFAULT_LOG_FORMAT = '[@%(hostname)s %(asctime)s] %(message)s'
coloredlogs.install(level='INFO')


class AnnotationPipeline:
    def __init__(self, cache, use_cache=True, use_web=True, proxies=None,
                 **cache_kwargs):
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
        - **cache_kwargs will be passed to the Cache constructor if the cache
          option is not already a Cache instance.

        See the docstring of Pipeline.run for some usage examples.
        """
        self.use_cache = use_cache
        self.use_web = use_web
        self.proxies = proxies

        if isinstance(cache, Cache):
            self.cache = cache
        else:
            self.cache = create_cache(cache, **cache_kwargs)

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
        msg = 'Starting annotation pipeline with options:\n\n{}\n'.format(opts)
        logger.info(msg)

        self._read_vcf(vcf_path)
        self._annotate_rs_variants()
        self._extract_entrez_gene_ids()
        self._annotate_genes()
        self._add_omim_variants_info()
        self._add_pubmed_entries_to_omim_entries()
        self._add_uniprot_variants_info()

        self.end_time = time.time()
        elapsed = format_timespan(self.end_time - self.start_time)
        logger.info('Done! Took {} to complete the pipeline'.format(elapsed))

        return self.rs_variants

    def _read_vcf(self, vcf_path):
        """Read the VCF and keep the variants with rs ID."""
        logger.info('Read "{}"'.format(vcf_path))
        self.variants = vcf_to_dataframe(vcf_path)
        have_single_rs = self.variants['id'].str.match(r'^rs\d+$')
        self.rs_variants = self.variants[have_single_rs].reset_index(drop=True)
        self.other_variants = self.variants[~have_single_rs].reset_index(drop=True)
        logger.info('{} variants with single rs'.format(len(self.rs_variants)))
        logger.info('{} other variants'.format(len(self.other_variants)))

    def _annotate_rs_variants(self):
        snp_annotator_classes = [
            ClinvarRsAnnotator,
            SnpeffAnnotator,
            MafAnnotator,
            HgvsAnnotator,
            DbsnpMyvariantAnnotator,
            DbsnpEntrezAnnotator
        ]

        for annotator_class in snp_annotator_classes:
            annotator = annotator_class(self.cache)
            annotator.PROXIES = self.proxies
            self._annotate_ids_and_add_series(
                    annotator=annotator,
                    id_series=self.rs_variants['id'],
                    df_to_modify=self.rs_variants,
                )

    def _extract_entrez_gene_ids(self):
        """Extract Entrez Gene IDs from dbsnp annotations. Adds a field in
        the dataframe self.rs_variants; returns a series of unique IDs."""

        def func(dbsnp_entries):
            if not dbsnp_entries: return []
            ids = {gene.get('geneid') for entry in dbsnp_entries
                                      for gene in entry.get('gene', {})}
            return list(ids)

        self.rs_variants['entrez_gene_ids'] = \
                self.rs_variants['dbsnp_myvariant'].fillna(False).map(func)

    def _annotate_genes(self):
        """Annotate Entrez genes; generate a dataframe with annotations."""
        # FIXME: gene annotation could follow the same logic than
        # omim variants and uniprot variants annotation later:
        # grab the unique ids and annotate in the same method
        # without adding an extra column in the rs_variants dataframe
        gene_ids = chain.from_iterable(self.rs_variants['entrez_gene_ids'])
        gene_annotations = pd.DataFrame({'entrez_id': list(set(gene_ids))})

        gene_annotator_classes = [MygeneAnnotator, GeneEntrezAnnotator]
        for annotator_class in gene_annotator_classes:
            annotator = annotator_class(cache=self.cache)
            self._annotate_ids_and_add_series(
                    annotator=annotator,
                    id_series=gene_annotations['entrez_id'],
                    df_to_modify=gene_annotations,
                )

        self.gene_annotations = gene_annotations

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
        for gene_variants in omim_variants_per_gene.values():
            for variant in gene_variants:
                rs = variant['rsid']
                if rs in self.rs_variants['id'].values:
                    # There might be more than one OMIM variant for a given rs
                    # --the typical case is the one of multiallelic SNPs.
                    # To simplify further parsing, we return the OMIM
                    # annotation for an rs always as a list:
                    if rs not in rs_to_omim_variants:
                        rs_to_omim_variants[rs] = []
                    rs_to_omim_variants[rs].append(variant)

        self.rs_variants['omim_entries'] = \
                self.rs_variants['id'].map(rs_to_omim_variants)

    def _add_pubmed_entries_to_omim_entries(self):
        pubmed_ids = set()

        for omim_entries in self.rs_variants['omim_entries'].dropna():
            for entry in omim_entries:
                for pubmed in entry.get('pubmeds', []):
                    pubmed_ids.add(pubmed['pmid'])

        pubmed_annotations = \
                PubmedAnnotator(cache=self.cache).annotate(pubmed_ids)

        for omim_entries in self.rs_variants['omim_entries'].dropna():
            for entry in omim_entries:
                for pubmed in entry.get('pubmeds', []):
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
        uniprot_variants = \
                pd.DataFrame(uniprot_variants).set_index('rsid', drop=False)

        def rs_to_uniprot_variants(rs):
            if rs not in uniprot_variants.index:
                return None
            return uniprot_variants.loc[[rs]].to_dict('records')

        self.rs_variants['uniprot_entries'] = \
            self.rs_variants['id'].map(rs_to_uniprot_variants)

    def _annotate_ids_and_add_series(self, annotator, id_series, df_to_modify):
        """
        Annotate the IDs in id_series with the passed annotator object and add
        the resulting Series to the df_to_modify.
        """
        annotations_dict = annotator.annotate(
                id_series,
                use_web=self.use_web,
                use_cache=self.use_cache,
            )
        annotations = id_series.map(annotations_dict)

        total = len(annotations)
        annotated_count = len(annotations[annotations.notnull()])
        msg = '{}/{} ({:.2%}) variants have {} data'
        logger.info(msg.format(annotated_count, total, annotated_count/total,
                               annotator.SOURCE_NAME))

        df_to_modify[annotator.SOURCE_NAME] = annotations

