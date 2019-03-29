from anotala.annotators.base_classes import MyVariantAnnotator
from anotala.helpers import infer_annotated_allele


class HgvsAnnotator(MyVariantAnnotator):
    SOURCE_NAME = 'hgvs'
    SCOPES = 'dbsnp.rsid'
    FIELDS = 'clinvar.hgvs,snpeff.ann.feature_id,snpeff.ann.hgvs_c,' + \
             'snpeff.ann.hgvs_p,emv.egl_protein,emv.egl_variant,evs.hgvs'

    @classmethod
    def _parse_hit(cls, hit):
        annotation = {'myvariant_hgvs_g': hit['_id']}
        annotation['genomic_allele'] = infer_annotated_allele(hit['_id'])

        annotation.update(cls._parse_hgvs_from_clinvar(hit))
        annotation.update(cls._parse_hgvs_from_snpeff(hit))
        annotation.update(cls._parse_hgvs_from_evs(hit))
        annotation.update(cls._parse_hgvs_from_emv(hit))

        return annotation

    @staticmethod
    def _parse_hgvs_from_evs(hit):
        """Given a myvariant annotation, try to extract HGVS from Exome
        Variant Service data. Returns a dict."""
        info = {}
        if 'evs' not in hit or 'hgvs' not in hit['evs']:
            return info

        hgvs = hit['evs']['hgvs']

        if 'coding' in hgvs:
            info['evs_hgvs_c'] = hgvs['coding']

        if 'protein' in hgvs:
            info['evs_hgvs_p'] = hgvs['protein']

        return info

    @staticmethod
    def _parse_hgvs_from_emv(hit):
        """Given a myvariant annotation, try to extract HGVS from EGL data in
        'EMV' key, whatever that is. Returns a dict."""
        info = {}
        if 'emv' not in hit:
            return info

        if 'egl_protein' in hit['emv']:
            info['egl_hgvs_p'] = hit['emv']['egl_protein']

        if 'egl_variant' in hit['emv']:
            info['egl_hgvs_c'] = hit['emv']['egl_variant']

        return info

    @staticmethod
    def _parse_hgvs_from_clinvar(hit):
        """
        Given a myvariant annotation, try to extract HGVS info from ClinVar
        data. Returns a dict.
        """
        info = {}
        if 'clinvar' not in hit:
            return info

        hgvs = hit['clinvar']['hgvs']

        if 'genomic' in hgvs:
            info['clinvar_hgvs_g'] = hgvs['genomic']

        if 'coding' in hgvs:
            info['clinvar_hgvs_c'] = hgvs['coding']

        if 'protein' in hgvs:
            info['clinvar_hgvs_p'] = hgvs['protein']

        return info

    @staticmethod
    def _parse_hgvs_from_snpeff(hit):
        """
        Given a myvariant annotation, try to extract HGVS info from SnpEff
        data. Returns a dict.
        """
        info = {}
        if 'snpeff' not in hit:
            return info

        info['snpeff_hgvs_c'] = []
        info['snpeff_hgvs_p'] = []

        annotations = hit['snpeff']['ann']
        # 'ann' is sometimes a list of annotations (each one a dict) and
        # sometimes a single dict. Make sure we always deal with lists:
        if isinstance(annotations, dict):
            annotations = [annotations]

        for annotation in annotations:
            feature = annotation.get('feature_id', '?')

            if 'hgvs_c' in annotation:
                hgvs_c = '{}:{}'.format(feature, annotation['hgvs_c'])
                info['snpeff_hgvs_c'].append(hgvs_c)

            if 'hgvs_p' in annotation:
                hgvs_p = '{}:{}'.format(feature, annotation['hgvs_p'])
                info['snpeff_hgvs_p'].append(hgvs_p)

        return info

