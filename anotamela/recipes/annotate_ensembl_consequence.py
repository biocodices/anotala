from more_itertools import always_iterable

from anotamela import EnsemblAnnotator
from anotamela.pipeline import extract_ensembl_consequence


def annotate_ensembl_consequence(rsids, cache='mysql'):
    """
    Given one or more rs IDs, annotate their ensembl 'most_severe_consequence'.
    Returns a dictionary of rsID => consequence.

    Accepts as *cache* either a cache name or a Cache object to initialize the
    EnsemblAnnotator with it.
    """
    ensembl = EnsemblAnnotator(cache=cache)
    annotations = ensembl.annotate(always_iterable(rsids))

    return {id_: extract_ensembl_consequence(annotation)
            for id_, annotation in annotations.items()}

