from anotamela.annotators.base_classes import AnnotatorWithCache


class LocalAnnotator(AnnotatorWithCache):
    """
    Base class for annotators that use a local database with annotations,
    instead of an online service. Classes that inherit from this one
    should have:

    - a class variable SOURCE_NAME, table of the database with the annotations.
      I intend to use this as the panel name, e.g. 'Nutripanel'.

    The class assumes the table will have these fields:

        - A "name" (typically an rs ID).
        - A "phenotype", JSON list of strings.
        - A "gene".
        - "alleles" (e.g. [A, G])
        - A "risk allele" (e.g. "A").
        - "full_annotations" with the full blown annotations for the risk
        - allele of that variant. These are meant for the investigator.
        - "report_annotations", with a selected annotation text that will be
          put in the reports.
        - A JSON dict field of "extra_info" that might have a "bibliography"
          key specific to that variant.
    """


