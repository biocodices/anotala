from anotamela.annotators.base_classes import MyVariantAnnotator


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
            # Make 'effects' always a list, even with only one prediction:
            effects = annotation.get('effect', '').split('&')
            annotation['effects'] = effects
            del(annotation['effect'])

        return annotations

