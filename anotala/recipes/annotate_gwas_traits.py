from more_itertools import always_iterable


from anotala import GwasCatalogAnnotator
from anotala.pipeline import extract_gwas_traits

def annotate_gwas_traits(rsids, **annotator_options):
    """
    Given one or more rs IDs, list the traits associated to them in GWAS.
    Returns a dictionary like { rs123: [trait-1, trait-2, ...] }

    Needs a *cache* that can be either a cache name ('myqsl') or a Cache
    instance, to initialize the annotator with it.
    """
    gwas = GwasCatalogAnnotator(**annotator_options)
    annotations = gwas.annotate(always_iterable(rsids))

    if annotations:
        return {id_: extract_gwas_traits(annotation)
                for id_, annotation in annotations.items()}
