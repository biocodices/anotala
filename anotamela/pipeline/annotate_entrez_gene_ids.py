from anotamela.annotators import ENTREZ_GENE_ANNOTATOR_CLASSES
from anotamela.pipeline import annotate_ids


def annotate_entrez_gene_ids(entrez_gene_ids, cache, use_web=True,
                             use_cache=True, proxies=None, sleep_time=None):
    """
    Given a list of *entrez_gene_ids*, return a dataframe with the passed IDs
    annotated with various annotators. The passed cache will be used for all
    the annotators: either a Cache subclass instance or a cache name
    ('postgres', 'redis') for each annotator to initialize.

    Extra kwargs are for the annotators:

        - use_cache (boolean) whether to use cached data
        - use_web (boolean) whether to use web annotation
        - proxies (dict) like {'http': 'socks5://localhost:9050'}
        - if not None, sleep_time will override each annotator's SLEEP_TIME
          as the time to sleep between queries.

    """
    annotator_classes = ENTREZ_GENE_ANNOTATOR_CLASSES.values()

    df = annotate_ids(entrez_gene_ids, annotator_classes, cache,
                      use_web=use_web, use_cache=use_cache, proxies=proxies,
                      sleep_time=sleep_time)
    df.rename(columns={'id': 'entrez_gene_id'}, inplace=True)
    return df

