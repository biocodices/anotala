from myvariant import MyVariantInfo

from anotamela.annotators.base_classes import AnnotatorWithCache


class MyVariantAnnotator(AnnotatorWithCache):
    """
    This class is meant as a class for annotators that use myvariant.info.
    To use it, define a new class that inherits from this one and that has:

        - SOURCE_NAME class variable
        - SCOPES class variable (scopes to query with the passed IDs)
        - FIELDS class variable (fields to retrieve for each ID)
        - an optional @staticmethod or @classmethod
          _parse_annotation(raw_annotation)
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
        for hit in hits:
            if 'notfound' not in hit:
                yield hit['query'], hit

