from anotamela.pipeline import annotate_ids
from anotamela.annotators import RS_ANNOTATOR_CLASSES


def annotate_rsids(rsids, cache, use_web=True, use_cache=True,
                   annotator_names='all', proxies={}, sleep_time=None):
    """
    Given a list of annotator names like 'clinvar', 'dbsnp_myvariant',
    'gwas_catalog', etc., and a list of *rsids*, return a dataframe with
    the passed IDs annotated with each of the passed annotators. If
    *annotators_names* is left as 'all', all the available annotators will be
    used.

    The passed cache will be used for all the annotators: either a Cache
    subclass instance or a cache name ('postgres', 'redis') for each annotator
    to initialize.

    The whole list of annotators available for rs IDs is found in:
    anotamela.annotators.RS_ANNOTATOR_CLASSES

    Extra kwargs are for the annotators:

        - use_cache (boolean) whether to use cached data
        - use_web (boolean) whether to use web annotation
        - proxies (dict) like {'http': 'socks5://localhost:9050'}
        - if not None, sleep_time will override each annotator's SLEEP_TIME
          as the time to sleep between queries.
    """
    if annotator_names == 'all':
        annotator_names = RS_ANNOTATOR_CLASSES.keys()

    annotator_classes = [RS_ANNOTATOR_CLASSES[name]
                         for name in annotator_names]

    df = annotate_ids(rsids, annotator_classes, cache, use_web=use_web,
                      use_cache=use_cache, proxies=proxies,
                      sleep_time=sleep_time)

    df.rename(columns={'id': 'rsid'}, inplace=True)
    return df
