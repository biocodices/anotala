from anotamela.annotators.base_classes import MyVariantAnnotator
from anotamela.helpers import infer_annotated_allele


class HgvsAnnotator(MyVariantAnnotator):
    SOURCE_NAME = 'hgvs'
    SCOPES = 'dbsnp.rsid'
    FIELDS = ('clinvar.hgvs snpeff.ann.feature_id snpeff.ann.hgvs_c '
              'snpeff.ann.hgvs_p').split()

    @staticmethod
    def _parse_hit(hit):
        annotation = {'myvariant_hgvs_g': hit['_id']}
        annotation['genomic_allele'] = infer_annotated_allele(hit['_id'])

        if 'snpeff' in hit:
            annotation['snpeff_hgvs_c'] = []
            annotation['snpeff_hgvs_p'] = []
            if isinstance(hit['snpeff']['ann'], dict):
                hit['snpeff']['ann'] = [hit['snpeff']['ann']]
            for ann in hit['snpeff']['ann']:
                feature = ann.get('feature_id', '?')
                if 'hgvs_c' in ann:
                    hgvs_c = '{}:{}'.format(feature, ann['hgvs_c'])
                    annotation['snpeff_hgvs_c'].append(hgvs_c)
                if 'hgvs_p' in ann:
                    hgvs_p = '{}:{}'.format(feature, ann['hgvs_p'])
                    annotation['snpeff_hgvs_p'].append(hgvs_p)

        if 'clinvar' in hit:
            annotation['clinvar_hgvs_g'] = hit['clinvar']['hgvs']['genomic']
            annotation['clinvar_hgvs_c'] = hit['clinvar']['hgvs']['coding']

        return annotation

