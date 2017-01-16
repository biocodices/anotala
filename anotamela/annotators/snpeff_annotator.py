from anotamela.annotators.base_classes import MyVariantAnnotator
from anotamela.helpers import (
    listify,
    infer_annotated_allele,
)


class SnpeffAnnotator(MyVariantAnnotator):
    SOURCE_NAME = 'snpeff_myvariant'
    SCOPES = 'dbsnp.rsid'
    FIELDS = 'snpeff'

    @classmethod
    def _parse_hit(cls, hit):
        annotations = hit['snpeff']['ann']

        # Make sure the annotations are always a *list* of annotations
        # Right now, myvariant client sometimes returns a list of annotaions
        # dictionaries and sometimes a single annotation dictionary.
        annotations = listify(annotations)

        for annotation in annotations:
            annotation['genomic_allele'] = hit['allele']
            annotation['coding_allele'] = \
                infer_annotated_allele(annotation.get('hgvs_c'))
            annotation['effects'] = annotation.get('effect', '').split('&')

            if 'effect' in annotation:
                del(annotation['effect'])

        return annotations

