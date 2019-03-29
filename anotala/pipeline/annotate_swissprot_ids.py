from itertools import chain

from anotala.annotators import UniprotAnnotator
from anotala.pipeline import annotate_ids


def annotate_swissprot_ids(swissprot_ids, cache, use_web=True, use_cache=True,
                           proxies={}, sleep_time=None):
    """
    Given a list of Swissprot IDs and annoator options, annotate them and
    return a flat list of all the variants found in those proteins.
    """
    df = annotate_ids(swissprot_ids, [UniprotAnnotator], cache=cache,
                      use_web=use_web, use_cache=use_cache, proxies=proxies,
                      sleep_time=sleep_time)

    return list(chain.from_iterable(df['uniprot'].dropna().values))

