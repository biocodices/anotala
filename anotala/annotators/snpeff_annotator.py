from anotala.annotators.base_classes import MyVariantAnnotator
from anotala.helpers import (
    listify,
    infer_annotated_allele,
    parse_prot_change,
)


class SnpeffAnnotator(MyVariantAnnotator):
    SOURCE_NAME = 'snpeff_myvariant'
    SCOPES = 'dbsnp.rsid'
    FIELDS = 'snpeff'

    @classmethod
    def _parse_hit(cls, hit):
        # Sometimes (e.g. rs203319) a SNP will have two alleles, one with
        # a SnpEff annotation and other without one, so this key will be
        # missing:
        if 'snpeff' not in hit:
            return None

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

            if 'hgvs_p' in annotation:
                annotation['hgvs_p'] = parse_prot_change(annotation['hgvs_p'])

        return annotations
