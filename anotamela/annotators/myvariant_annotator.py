from myvariant import MyVariantInfo


class MyVariantAnnotator():
    """This class is meant as a second parent class for AnnotatorWithCache
    subclasses that use myvariant.info to get the data."""

    def annotate(self, ids, use_cache=True, use_web=True, parse_data=True):
        # This wrapper around AnnotatorWithCache.annotate is meant to remove
        # the parallel and sleep_time arguments from annotators that inherit
        # from MyVariantAnnotator. Those options are superflous for them since
        # the myvariant client already deals with batch queries.
        return super().annotate(ids, parallel=len(ids), sleep_time=0,
                                use_cache=use_cache, use_web=use_web,
                                parse_data=parse_data)

    def _myvariant_query_and_cache(self, ids, scopes, fields):
        """Uses myvariant.info service to query many IDs across the given
        scopes. It returns a dict of {id: result, ... } with the IDs that
        were found (i.e. leaves out the not found ones). The successful
        results are also cached."""
        if not hasattr(self, 'mv'):
            self.mv = MyVariantInfo()
        hits = self.mv.querymany(ids, scopes=scopes, fields=fields)
        annotations = {hit['query']: hit for hit in hits
                       if 'notfound' not in hit}
        self.cache.set(annotations, namespace=self.SOURCE_NAME)
        return annotations

