import re

from anotamela.annotators.base_classes import MyVariantAnnotator
from anotamela.helpers import (
    listify,
    infer_annotated_allele,
    parse_prot_change,
)


def _build_variant_regex():
    """
    Helper method for ClinvarRsAnnotator: it compiles regular expressions
    that will be used to extract the aminoacid change, transcript name,
    cds change and gene symbol from the 'preferred_name' that ClinVar gives
    to a variant.
    """
    patterns = {
        'cds': r'(?P<cds_change>c\..+)',
        'gene': r'(?P<gene>.+)',
        'prot_change': r'\((?P<prot_change>p\..+)\)',
        'transcript': r'(?P<transcript>NM_.+)',
    }
    patterns['gene_in_parens'] = r'\({gene}\)'.format(**patterns)
    patterns['full_name'] = \
        r'{transcript}{gene_in_parens}:{cds} {prot_change}'.format(**patterns)
    patterns['no_protchange'] = \
        r'{transcript}{gene_in_parens}:{cds}'.format(**patterns)
    patterns['no_transcript'] = \
        r'{gene}:{cds} {prot_change}'.format(**patterns)

    pattern_names = 'full_name no_protchange no_transcript'.split()
    return {pattern_name: re.compile(patterns[pattern_name])
            for pattern_name in pattern_names}


class ClinvarRsAnnotator(MyVariantAnnotator):
    SOURCE_NAME = 'clinvar_rs'
    SCOPES = 'clinvar.rsid'
    FIELDS = 'clinvar'
    VARIANT_REGEX = _build_variant_regex()

    @classmethod
    def _parse_hit(cls, hit):
        if 'clinvar' not in hit:
            return []

        rcv_annotations = listify(hit['clinvar']['rcv'])

        for annotation in rcv_annotations:
            # Some significances are 'compound', like "Pathogenic, risk factor"
            # We parse them always to a list like ['Pathogenic', 'risk factor']:
            annotation['clinical_significances'] = \
                annotation.get('clinical_significance', '').split(', ')

            if 'clinical_significance' in annotation:
                del(annotation['clinical_significance'])

            annotation['conditions'] = listify(annotation['conditions'])

            for condition in annotation['conditions']:
                condition.update(cls._parse_condition_identifiers(condition))
                if 'identifiers' in condition:
                    del(condition['identifiers'])

            annotation['url'] = cls._url(annotation['accession'])
            annotation['variant_url'] = cls._variant_url(hit['clinvar']['variant_id'])

            name_info = cls._parse_preferred_name(annotation['preferred_name'])
            annotation.update(name_info)

            if 'prot_change' in annotation:
                annotation['prot_change'] = \
                    parse_prot_change(annotation['prot_change'])

            annotation['genomic_allele'] = hit['allele']
            annotation['coding_allele'] = \
                infer_annotated_allele(annotation.get('cds_change', ''))

            # Copy the variant info to each particular RCV entry:
            for key, value in hit['clinvar'].items():
                if key == 'rcv':
                    continue

                # Flatten the dicts one level
                if isinstance(value, dict):
                    for key2, val2 in value.items():
                        compound_key = '{}_{}'.format(key, key2)
                        annotation[compound_key] = val2
                    continue

                annotation[key] = value

        return rcv_annotations

    @staticmethod
    def _url(accession):
        return f'https://www.ncbi.nlm.nih.gov/clinvar/{accession}'

    @staticmethod
    def _variant_url(variant_id):
        return f'https://www.ncbi.nlm.nih.gov/clinvar/variation/{variant_id}'

    @classmethod
    def _parse_preferred_name(cls, preferred_name):
        matches_1 = cls.VARIANT_REGEX['full_name'].search(preferred_name)
        matches_2 = cls.VARIANT_REGEX['no_protchange'].search(preferred_name)
        matches_3 = cls.VARIANT_REGEX['no_transcript'].search(preferred_name)

        if matches_1:
            return matches_1.groupdict()
        if matches_2:
            return matches_2.groupdict()
        if matches_3:
            return matches_3.groupdict()

        return {}

    @staticmethod
    def _parse_condition_identifiers(condition):
        """
        Given a ClinVar condition dictionary, relocate IDs found under
        'identifiers' to the root level, with a key that indicates their
        origin. E.g. {'identifiers': {'omim': 1}} becomes {'omim_id': 1}
        """
        ids = {}

        if 'identifiers' not in condition:
            return ids

        for resource_name, id_ in condition['identifiers'].items():
            new_key = '{}_id'.format(resource_name)
            ids[new_key] = id_

        return ids

