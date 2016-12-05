import pandas as pd

from anotamela.annotators import (AnnotatorWithCache,
                                  OmimGeneAnnotator)


class OmimVariantAnnotator(AnnotatorWithCache):
    """
    Provider of OMIM annotations for gene variants (e.g. 605557.0001).
    Responses are cached.

        > omim_variant_annotator = OmimVariantAnnotator()
        > omim_variant_annotator.annotate('605557.0001')
        # => { '605557.0001': ... }

    The base sleep time between requests is 60s because OMIM is quite strict
    against crawlers. However, you can lower this value at your own risk:

        > omim_variant_annotator.SLEEP_TIME = 30

    SLEEP_TIME is used as a base value that will be randomized each time.
    The user agent of each request is also randomized.
    """
    SOURCE_NAME = 'omim_variants'
    ANNOTATIONS_ARE_JSON = True

    SLEEP_TIME = 60

    def _batch_query_and_cache(self, ids, _, __):
        # This annotator uses OmimGeneAnnotator to get all the variants in a
        # given gene, and then extract the variants that were requested.
        # The gene annotation step is necessary because the variant info is in
        # the gene entry page, which has to be parsed completely before we can
        # extract the variant bit. In the process, the OmimGeneAnnotator
        # caches its own gene raw data from OMIM.
        if not hasattr(self, 'omim_gene_annotator'):
            self.omim_gene_annotator = OmimGeneAnnotator(cache=self.cache)

        gene_ids = [id_.split('.')[0] for id_ in ids]
        gene_dataframes = self.omim_gene_annotator.annotate(gene_ids)
        df = pd.concat(gene_dataframes.values()).set_index('variant_id')
        # Check each variant has a single row in the dataframe:
        assert len(df.index) == len(set(df.index))

        # From all the gene variants, only keep the ones that were queried
        annotations = df.loc[ids].to_dict('index')
        self.cache.set(annotations, namespace=self.SOURCE_NAME)
        return annotations

