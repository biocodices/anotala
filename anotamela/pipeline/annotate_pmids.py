from anotamela.annotators import PubmedAnnotator
from anotamela.pipeline import annotate_ids


def annotate_pmids(pmids, cache, use_web=True, use_cache=True, proxies={},
                   sleep_time=None):
    """
    Given a list of PubMed IDs and annotator options, annotate them and return
    a dictionary of annotations by pmid.
    """
    df = annotate_ids(pmids, [PubmedAnnotator], cache=cache, use_web=use_web,
                      use_cache=use_cache, proxies=proxies,
                      sleep_time=sleep_time)

    df.rename(columns={'id': 'pmid'}, inplace=True)
    return df.set_index('pmid', drop=False).to_dict('index')

