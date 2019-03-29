from more_itertools import always_iterable

from anotala import EnsemblAnnotator
from anotala.pipeline import extract_ensembl_consequence


def annotate_ensembl_consequence(rsids, **annotator_options):
    """
    Given one or more rs IDs, annotate their ensembl 'most_severe_consequence'.
    Returns a dictionary of rsID => consequence.

    Accepts a *cache* that can be either a cache name ('myqsl') or a Cache
    instance, to initialize the annotator with it.
    """
    ensembl = EnsemblAnnotator(**annotator_options)
    annotations = ensembl.annotate(always_iterable(rsids))

    if annotations:
        return {id_: extract_ensembl_consequence(annotation)
                for id_, annotation in annotations.items()}

