from anotamela.annotators import AnnotatorWithCache, MyVariantAnnotator


class HgvsAnnotator(MyVariantAnnotator, AnnotatorWithCache):
    SOURCE_NAME = 'hgvs'

    def _batch_query(self, ids, _, __):
        hgvs_fields = ('clinvar.hgvs snpeff.ann.feature_id snpeff.ann.hgvs_c '
                       'snpeff.ann.hgvs_p').split()
        return self._myvariant_query_and_cache(ids, scopes='dbsnp.rsid',
                                               fields=hgvs_fields)

    @staticmethod
    def _parse_annotation(raw_annotation):
        d = raw_annotation
        annotation = {'myvariant_hgvs_g': d['_id']}

        if 'snpeff' in d:
            annotation['snpeff_hgvs_c'] = []
            annotation['snpeff_hgvs_p'] = []
            if isinstance(d['snpeff']['ann'], dict):
                d['snpeff']['ann'] = [d['snpeff']['ann']]
            for ann in d['snpeff']['ann']:
                feature = ann.get('feature_id', '?')
                if 'hgvs_c' in ann:
                    hgvs_c = '{}:{}'.format(feature, ann['hgvs_c'])
                    annotation['snpeff_hgvs_c'].append(hgvs_c)
                if 'hgvs_p' in ann:
                    hgvs_p = '{}:{}'.format(feature, ann['hgvs_p'])
                    annotation['snpeff_hgvs_p'].append(hgvs_p)

        if 'clinvar' in d:
            annotation['clinvar_hgvs_g'] = d['clinvar']['hgvs']['genomic']
            annotation['clinvar_hgvs_c'] = d['clinvar']['hgvs']['coding']

        return annotation

