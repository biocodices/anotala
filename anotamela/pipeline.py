import time
import logging
from itertools import chain

import pandas as pd
import coloredlogs
from vcf_to_dataframe import vcf_to_dataframe

from anotamela.cache import create_cache
from anotamela.annotators import *


logger = logging.getLogger(__name__)
coloredlogs.DEFAULT_LOG_FORMAT = '[@%(hostname)s %(asctime)s] %(message)s'
coloredlogs.install(level='INFO')


class Pipeline:
    def run(self, vcf_path, out_filepath=None, use_cache=True, use_web=True,
            proxies=None, cache=None, **cache_kwargs):
        """
        Annotate the given VCF file (accepts gzipped file too). Return the
        filepath to a CSV with annotations per variant. Options are:

        - out_filepath (default=None): path where the output CSV with
          annotations should be written. If None, it will use the same filepath
          as the input VCF replacing '.vcf[.gz]' with '.csv'.
        - use_cache (default=True): whether to use or not cached data about
          variant.
        - use_web (default=True): whether to use or not web data to annotate
          the variants. If use_cache is also set, the web will be used only
          to annotate the ones not found in cache.
        - proxies (default=None) can be a dictionary of proxies that will be
          used by the requests library for some annotators. For instance:
          {'http': 'socks5://localhost:9050'}
        - cache (default=None): you can 'manually' instantiate a Cache
          (either PostgresCache or RedisCache and pass it here), or you can
          specify 'redis' or 'postgres' and let the pipeline do that for you.
        - **cache_kwargs will be passed to the Cache constructor if the cache
          option is not already a Cache instance.

        Examples:

            > from anotamela import run_pipeline
            > run_pipeline('~/variants.vcf.gz',
                           out_filepath='~/annotations.csv')

            > from anotamela import run_pipeline
            > run_pipeline('~/variants.vcf.gz', cache='redis',
                           out_filepath='~/annotations.csv')

            > from anotamela import run_pipeline
            > run_pipeline('~/variants.vcf', cache='postgres',
                           credentials_filepath='~/.pg_creds.yml')

            > from anotamela import run_pipeline
            > run_pipeline('~/variants.vcf', cache='redis',
                           use_cache=False, host='192.168.1.10')
            # Will reannotate any already cached variants, but it will cache
            # the new annotations, so the Redis connection will be used

        """
        self.start_time = time.time()

        # Options
        self.cache = create_cache(cache, **cache_kwargs)
        self.use_cache = use_cache
        self.use_web = use_web
        self.proxies = proxies

        # Pipeline
        self._read_vcf(vcf_path)
        self._annotate_snps(self.rs_variants)
        gene_ids = self._extract_gene_ids()
        self.gene_annotations = self._annotate_genes(gene_ids)
        self._add_omim_variants_info()
        self._add_pubmed_entries_to_omim_entries()
        self._add_uniprot_variants_info()

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

    def _add_omim_variants_info(self):
        annotator = OmimGeneAnnotator(cache=self.cache)
        annotator.PROXIES = self.proxies
        frames = annotator.annotate_from_entrez_ids(
                self.gene_annotations['entrez_id'],
                use_cache=self.use_cache,
                use_web=self.use_web,
            )
        omim_variants = pd.concat(frames, ignore_index=True)
        omim_variants = omim_variants.set_index('rsid', drop=False)\
                                     .loc[self.rs_variants['id']]
        omim_variants = omim_variants.dropna(how='all')

        def rs_to_omim_variants(rs):
            if rs not in omim_variants.index:
                return None
            # .loc[[rs]] forces the result to always be a DataFrame, whereas
            # .loc[rs] will give you a Series if there's one row, or a Df if
            # there are many matching rows:
            return omim_variants.loc[[rs]].to_dict('records')

        self.rs_variants['omim_entries'] = \
                self.rs_variants['id'].map(rs_to_omim_variants)

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
        uniprot_variants = annotator.annotate(
            uniprot_ids,
            use_cache=self.use_cache,
            use_web=self.use_seb
        )
        uniprot_variants = uniprot_variants.set_index('rsid', drop=False)

        def rs_to_uniprot_variants(rs):
            if rs not in uniprot_variants.index:
                return None
            return uniprot_variants.loc[[rs]].to_dict('records')

        self.rs_variants['uniprot_entries'] = \
            rs_variants['id'].map(rs_to_uniprot_variants)

        return out_filepath

    def _annotate_genes(self, gene_ids):
        """Annotate Entrez genes; generate a dataframe with annotations."""
        gene_annotations = pd.DataFrame({'entrez_id': gene_ids})

        gene_annotator_classes = [MygeneAnnotator, GeneEntrezAnnotator]
        for annotator_class in gene_annotator_classes:
            annotator = annotator_class(cache=self.cache)
            self._annotate_ids_and_add_series(
                    annotator=annotator,
                    id_series=gene_annotations['entrez_id'],
                    df_to_modify=gene_annotations,
                )

        return gene_annotations

        self.gene_annotations = gene_annotations

    def _read_vcf(self, vcf_path):
        """Read the VCF and keep the variants with rs ID."""
        logger.info('Read "{}"'.format(vcf_path))
        self.variants = vcf_to_dataframe(vcf_path)
        have_single_rs = self.variants['id'].str.match(r'^rs\d+$')
        self.rs_variants = self.variants[have_single_rs].reset_index(drop=True)
        self.other_variants = self.variants[~have_single_rs].reset_index(drop=True)
        logger.info('{} variants with single rs'.format(len(self.rs_variants)))
        logger.info('{} other variants'.format(len(self.other_variants)))

    def _annotate_snps(self, rs_variants_df):
        snp_annotator_classes = [
            ClinvarRsAnnotator,
            #  SnpeffAnnotator,
            #  MafAnnotator,
            #  HgvsAnnotator,
            DbsnpMyvariantAnnotator,
            #  DbsnpEntrezAnnotator
        ]

        for annotator_class in snp_annotator_classes:
            annotator = annotator_class(self.cache)
            annotator.PROXIES = self.proxies
            self._annotate_ids_and_add_series(
                    annotator=annotator,
                    id_series=rs_variants_df['id'],
                    df_to_modify=rs_variants_df,
                )

    def _extract_gene_ids(self):
        """Extract Entrez Gene IDs from dbsnp annotations. Adds a field in
        the dataframe self.rs_variants; returns a series of unique IDs."""
        gene_ids_per_rs = self.rs_variants['dbsnp_myvariant'].fillna(False).map(
                self._extract_gene_id_from_dbsnp_entries
            )
        self.rs_variants['entrez_gene_ids'] = gene_ids_per_rs
        return list(set(chain.from_iterable(gene_ids_per_rs)))

    def _annotate_ids_and_add_series(self, annotator, id_series, df_to_modify):
        """
        Given a DataFrame with a column of IDs by the name of <id_series>,
        annotate those IDs with the passed annotator object (it will call
        annotator.annotate() with the passed **annotate_kwargs if any).

        The resulting annotations will be added as a new column to the
        dataframe, by the name of annotator.SOURCE_NAME
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

    @staticmethod
    def _extract_gene_id_from_dbsnp_entries(dbsnp_entries):
        if not dbsnp_entries:
            return []
        ids = {gene.get('geneid') for entry in dbsnp_entries
                                for gene in entry.get('gene', {})}
        return list(ids)

    @staticmethod
    def _mims_from_rcvs(rcvs):
        """Extract OMIM variant IDs from a RCV ClinVar list."""
        rcvs = rcvs or []
        mim_ids = {rcv['omim'] for rcv in rcvs if 'omim' in rcv}
        if not mim_ids:
            return
        assert len(mim_ids) == 1  # Should have 1 or 0 MIM IDs per rs variant
        return mim_ids.pop()

