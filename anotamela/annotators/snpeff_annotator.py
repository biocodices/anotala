from anotamela.annotators.base_classes import MyVariantAnnotator
from anotamela.helpers import SNP_RE


class SnpeffAnnotator(MyVariantAnnotator):
    SOURCE_NAME = 'snpeff_myvariant'
    SCOPES = 'dbsnp.rsid'
    FIELDS = 'snpeff'

    @staticmethod
    def _parse_hit(hit):
        annotations = hit['snpeff']['ann']

        # Make sure the annotations are always a *list* of annotations
        # Right now, myvariant client sometimes returns a list of annotaions
        # dictionaries and sometimes a single annotation dictionary.
        if not isinstance(annotations, list):
            annotations = [annotations]

        for annotation in annotations:
            # Make 'effects' always a list, even with one or zero predictions
            effects = annotation.get('effect', '').split('&')
            annotation['effects'] = effects
            del(annotation['effect'])

            # Infer the allele being describe in the annotation
            snp_match = SNP_RE.search(annotation['hgvs_c'])
            if snp_match:
                annotation['allele'] == snp_match.group('new_allele')
            else:
                annotation['allele'] == annotation['hgvs_c']

        return annotations

