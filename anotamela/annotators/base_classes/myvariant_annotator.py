from operator import itemgetter
from itertools import groupby

from myvariant import MyVariantInfo

from anotamela.annotators.base_classes import AnnotatorWithCache


class MyVariantAnnotator(AnnotatorWithCache):
    """
    This class is meant as a class for annotators that use myvariant.info.
    To use it, define a new class that inherits from this one and that has:

        - SOURCE_NAME class variable
        - SCOPES class variable (scopes to query with the passed IDs)
        - FIELDS class variable (fields to retrieve for each ID)
        - a @staticmethod or @classmethod _parse_hit(hit) to parse each of the
          hits (results) associated to a single ID
        - optionally, set a class variable VERBOSE to get all the output
          produced by myvariant (silenced by default)
    """
    ANNOTATIONS_ARE_JSON = True
    VERBOSE = False

    def _batch_query(self, ids):
        """
        Uses myvariant.info service to query many IDs across the scopes defined
        by the class variable SCOPES. It returns a dict of {id: result, ... }
        with the IDs that were found (i.e. leaves out the not found ones) and
        info from the fields defined in the class variable FIELDS.

        The successful results are cached.
        """
        if not hasattr(self, 'mv'):
            self.mv = MyVariantInfo()

        hits = self.mv.querymany(ids, scopes=self.SCOPES, fields=self.FIELDS,
                                 verbose=self.VERBOSE)

        for query, hits_group in groupby(hits, itemgetter('query')):
            hits_group = list(hits_group)
            if 'notfound' not in hits_group[0]:
                # Yields a tuple (ID, group of annotations for the ID)
                yield query, hits_group

    @classmethod
    def _parse_annotation(cls, hits_group):
        annotations = [cls._parse_hit(hit) for hit in hits_group]
        if annotations:
            return annotations

    @staticmethod
    def _parse_hit(hit):
        raise NotImplementedError

