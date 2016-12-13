import re

from anotamela.annotators.base_classes import MyVariantAnnotator
from anotamela.helpers import listify


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
    def _parse_annotation(cls, hits):
        """
        Transforms each of the ClinVar annotations returned by MyVariant to a
        list of RCV entries. Each RCV describes an association of the given rs
        ID to a condition. The condition will vary between items in the list,
        but the variant info will be repeated (e.g. rsid, omim variant id,
        gene affected, genomic coordinates).

        Overrides MyvariantAnnotator._parse_annotation to add the flattening
        of the list of lists logic.
        """
        rcv_list = []
        for hit in hits:
            if 'clinvar' in hit:
                rcv_list += cls._parse_hit(hit)

        if rcv_list:
            return rcv_list

    @classmethod
    def _parse_hit(cls, hit):
        try:
            rcvs = listify(hit['clinvar']['rcv'])
        except TypeError:
            import q; q(hit)
            raise

        for rcv in rcvs:
            rcv['conditions'] = listify(rcv['conditions'])

            for condition in rcv['conditions']:
                if 'identifiers' in condition:
                    for resource_name, id_ in condition['identifiers'].items():
                        new_key = '{}_id'.format(resource_name)
                        condition[new_key] = id_
                    del(condition['identifiers'])

            rcv['url'] = cls.url(rcv['accession'])
            rcv.update(cls._parse_preferred_name(rcv['preferred_name']))

            # Copy the variant info to each particular RCV entry:
            for key, value in hit['clinvar'].items():
                if key == 'rcv':
                    continue
                # Flatten the dicts one level
                if isinstance(value, dict):
                    for key2, val2 in value.items():
                        compound_key = '{}_{}'.format(key, key2)
                        rcv[compound_key] = val2
                    continue

                rcv[key] = value

        return rcvs

    @staticmethod
    def url(accession):
        return 'https://www.ncbi.nlm.nih.gov/clinvar/{}'.format(accession)

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

