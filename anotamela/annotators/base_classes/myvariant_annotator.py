from myvariant import MyVariantInfo

from anotamela.annotators import AnnotatorWithCache


class MyVariantAnnotator(AnnotatorWithCache):
    """
    This class is meant as a class for AnnotatorWithCache
    subclasses that use myvariant.info to get the data.

    To use it, define a new class that inherits from this one and that has:

        - SOURCE_NAME class variable
        - SCOPES class variable (scopes to query)
        - FIELDS class variable (fields to retrieve)
        - an optional @staticmethod _parse_annotation(raw_annotation)
        - optionally, set <class>.VERBOSE = True to get all the output produced
          by myvariant.
    """
    ANNOTATIONS_ARE_JSON = True
    VERBOSE = False

    def _batch_query_and_cache(self, ids):
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
        annotations = {hit['query']: hit for hit in hits
                       if 'notfound' not in hit}
        self.cache.set(annotations, namespace=self.SOURCE_NAME)
        return annotations

