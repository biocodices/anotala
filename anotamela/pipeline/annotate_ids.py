import logging

import pandas as pd


logger = logging.getLogger(__name__)


def annotate_ids(ids, annotator_classes, cache, use_web=True, use_cache=True,
                 proxies={}, sleep_time=None):
    """
    Given a list of *ids* and *annotator_classes*, annotate all the ids with
    all the given annotators, using the passed *cache*. The rest of the
    kwargs apply to all the annotators.

    This method is meant as a common logic for annotate_rsids and
    annotate_entrez_ids methods.
    """
    df = pd.DataFrame({'id': ids})

    for annotator_class in annotator_classes:
        annotator = annotator_class(cache)

        if proxies:
            annotator.PROXIES = proxies

        if sleep_time is not None:
            annotator.SLEEP_TIME = sleep_time

        annotations = annotator.annotate(df['id'], use_web=use_web,
                                         use_cache=use_cache)
        df[annotator.SOURCE_NAME] = df['id'].map(annotations)

        found = len(annotations)
        total = len(df)
        logger.info('{}/{} ({:.2%}) variants have {} data'
                    .format(found, total, found/total, annotator.SOURCE_NAME))

    return df



