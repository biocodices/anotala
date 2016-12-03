from myvariant import MyVariantInfo

from anotamela.annotators import AnnotatorWithCache


class ClinvarRsAnnotator(AnnotatorWithCache):
    SOURCE_NAME = 'clinvar_rs'

    def annotate(self, ids, use_cache=True, use_web=True):
        return super().annotate(ids, parallel=len(ids), sleep_time=0,
                                use_cache=use_cache, use_web=use_web,
                                parse_data=False)

    def _batch_query(self, ids, _, __):
        mv = MyVariantInfo()
        hits = mv.querymany(ids, scopes='clinvar.rsid', fields='clinvar')
        return {hit['query']: hit['clinvar'] for hit in hits
                if 'clinvar' in hit}

