class MyVariantAnnotator():
    """This class is meant as a second parent class for AnnotatorWithCache
    subclasses that also happen to use myvariant to get the data."""

    def annotate(self, ids, use_cache=True, use_web=True, parse_data=True):
        # This wrapper around AnnotatorWithCache.annotate is meant to remove
        # the parallel and sleep_time arguments from annotators that inherit
        # from MyVariantAnnotator, because those options are superflous for
        # them since the myvariant client already deals with batch queries.
        return super().annotate(ids, parallel=len(ids), sleep_time=0,
                                use_cache=use_cache, use_web=use_web,
                                parse_data=parse_data)

