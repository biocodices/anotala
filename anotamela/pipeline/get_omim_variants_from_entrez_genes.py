from anotamela import OmimGeneAnnotator


def get_omim_variants_from_entrez_genes(entrez_gene_ids, cache,
                                        use_web=True, use_cache=True,
                                        proxies=None, sleep_time=None):
    """
    Given a list of Entrez genes, annotate the OMIM genes IDs associated to
    them and return a dictionary with the variants found in those genes and
    rs IDs as keys. Any variants without an rs ID are not included.

    Kwargs are for the annotators:
        - cache, either a cache name like 'postgres' or a Cache instance
        - use_web, use_cache, options for the annotation
        - proxies for the annotation via Omim scraping, typically some Tor
          instance, described as {'http': 'socks5://localhost:9050'}
        - sleep_time: optional time in seconds to sleep between queries. If not
          set, the default SLEEP_TIME of each OmimGeneAnnotator will be used.

    """
    # Annotate with OMIM
    annotator = OmimGeneAnnotator(cache=cache, proxies=proxies)
    if sleep_time:
        annotator.SLEEP_TIME = sleep_time
    omim_variants = annotator.annotate_from_entrez_ids(entrez_gene_ids,
                                                       use_cache=use_cache,
                                                       use_web=use_web)

    # Get all the variants in the same flat list:
    omim_variants = [entry for gene_variants in omim_variants.values()
                           for entry in gene_variants]

    return omim_variants

