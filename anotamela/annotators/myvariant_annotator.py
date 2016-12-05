from myvariant import MyVariantInfo


class MyVariantAnnotator():
    """
    This class is meant as a second parent class for AnnotatorWithCache
    subclasses that use myvariant.info to get the data.

    To use it, define a new class that inherits from this one and from
    AnnotatorWithCache and that has:
        - SOURCE_NAME class variable
        - SCOPES class variable (scopes to query)
        - FIELDS class variable (fields to retrieve)
        - an optional @staticmethod _parse_annotation(raw_annotation)
    """
    ANNOTATIONS_ARE_JSON = True

    def annotate(self, ids, use_cache=True, use_web=True, parse_data=True):
        # This wrapper around AnnotatorWithCache.annotate is meant to remove
        # the parallel and sleep_time arguments from annotators that inherit
        # from MyVariantAnnotator. Those options are superflous for them since
        # the myvariant client already deals with batch queries.
        return super().annotate(ids, parallel=len(ids), sleep_time=0,
                                use_cache=use_cache, use_web=use_web,
                                parse_data=parse_data)

    def _batch_query_and_cache(self, ids, _, __):
        # The ignored arguments _ and __ are there to handle the <parallel>
        # and <sleep_time> arguments that AnnotatorWithCache.annotate
        # assumes.
        # Uses myvariant.info service to query many IDs across the given
        # scopes. It returns a dict of {id: result, ... } with the IDs that
        # were found (i.e. leaves out the not found ones). The successful
        # results are also cached.
        if not hasattr(self, 'mv'):
            self.mv = MyVariantInfo()
        hits = self.mv.querymany(ids, scopes=self.SCOPES, fields=self.FIELDS)
        annotations = {hit['query']: hit for hit in hits
                       if 'notfound' not in hit}
        self.cache.set(annotations, namespace=self.SOURCE_NAME)
        return annotations

