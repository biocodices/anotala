from anotamela.annotators.base_classes import MyVariantAnnotator
from anotamela.helpers import (
    SNP_RE,
    DEL_RE,
    listify
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
            annotation['allele'] = cls._infer_allele(annotation)
            annotation['effects'] = cls._split_effects(annotation)

            if 'effect' in annotation:
                del(annotation['effect'])

        return annotations

    @staticmethod
    def _infer_allele(annotation):
        """Given a SnpEff annotation, infer the allele that is described."""
        allele = None

        if not 'hgvs_c' in annotation:
            return allele

        snp_match = SNP_RE.search(annotation['hgvs_c'])
        del_match = DEL_RE.search(annotation['hgvs_c'])

        if snp_match:
            allele = snp_match.group('new_allele')
        elif del_match:
            allele = del_match.group('new_allele')
        else:
            # If extractions fail, return the whole change
            allele = annotation['hgvs_c']

        return allele

    @staticmethod
    def _split_effects(annotation):
        """
        Given a SnpEff annotation, split the 'effects' described into a list
        or return an empty list if no effects are described.
        """
        effects = []

        if 'effect' in annotation:
            effects += annotation['effect'].split('&')

        return effects

