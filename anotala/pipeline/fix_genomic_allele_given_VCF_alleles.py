import re
from copy import deepcopy
import logging

import pandas as pd


logger = logging.getLogger(__name__)


def fix_genomic_alleles_for_variant(variant):
    """
    Given a variant (pd.Series or dict), fix the genomic alleles at all
    columns/keys where it's possible. The fix happens inplace!
    """
    keys_to_skip = [
        'chrom',
        'pos',
        'id',
        'ref',
        'alt',
        'qual',
        'filter',
        'info',
        'format',
        'position_tag',
        'entrez_gene_ids',
        'entrez_gene_symbols',
    ]
    for key in variant.index:
        if key in keys_to_skip:
            continue
        try:
            variant[key] = fix_genomic_allele_given_VCF_alleles(
                entry_or_entries=variant[key],
                ref=variant['ref'], alts=variant['alt']
            )
        except ValueError:
            logger.warning(f'Genomic allele fix failed at key: "{key}"')

    return variant


def fix_genomic_allele_given_VCF_alleles(entry_or_entries, ref, alts):
    """
    Given either a single entry or a list of entries, each a dictionary with
    a 'genomic_allele', return a copy of the entry_or_entries with the genomic
    allele "fixed" to match the alleles from the VCF (ref, alts). This is
    meant to solve the way that some annotators name indels vs. the way
    the VCF names indels, e.g. "insC" vs "AC" or "del" vs "A".

    Returns a single dict if the input was a dict, or a list of dicts if the
    input was a list of dicts.
    """
    if isinstance(entry_or_entries, list):
        fixed = [_fix_genomic_allele_for_single_entry(entry, ref, alts)
                 for entry in entry_or_entries]
    elif isinstance(entry_or_entries, dict):
        fixed = _fix_genomic_allele_for_single_entry(entry_or_entries, ref, alts)
    elif pd.isnull(entry_or_entries):
        fixed = None
    else:
        raise ValueError('I was expecting a dict or a list of dicts, ' +
                         f'but I got a {type(entry_or_entries)}: ' +
                         str(entry_or_entries))

    return fixed


def _fix_genomic_allele_for_single_entry(entry, ref, alts):
    """Auxiliary function.
    See fix_genomic_allele_given_VCF_alleles for the explanation."""
    if not isinstance(entry, dict):
        raise ValueError(f'I was expecting a dict, I got a {type(entry)}: ' +
                         str(entry))

    entry_allele = entry.get('genomic_allele')

    if not entry_allele:
        return entry

    alleles = set([ref] + alts)

    # Do not fix an allele that's already as the VCF alleles
    if entry_allele in alleles:
        return entry

    # VCF notation for indel includes the nucleotide previous to the mutation
    # itself. So, for an insertion of a "C" after an "A", the alleles are
    # "A" and "AC", and for a deletion of a "C" after an "A", the alleles are
    # "AC", and "A". We need to match something like "insC" or "del" to "AC"
    # and "A" respectively, detecting the nucleotide at -1 and adding it.
    is_indel = any(len(allele) > 1 for allele in alleles)
    common_previous_nucleotides = set(allele[0] for allele in alleles)

    # NOTE: multiallelic deletions are problematic to choose an allele:
    # from "del" given "ACCC" I can't tell if the allele is "ACC" vs "AC".
    multiallelic_del = ('del' in entry_allele) and len(alleles) > 2

    if is_indel and len(common_previous_nucleotides) != 1:
        logger.warning(f"{alleles} don't have a unique previous nucleotide?")
        common_previous_nucleotide = None
    else:
        common_previous_nucleotide = common_previous_nucleotides.pop()

    fixed_entry = deepcopy(entry)
    if entry_allele and is_indel and common_previous_nucleotide and not multiallelic_del:
        fixed_allele = re.sub(r'ins|del', '', entry_allele)
        fixed_allele = f'{common_previous_nucleotide}{fixed_allele}'

        if fixed_allele in alleles:
            fixed_entry['genomic_allele'] = fixed_allele

            # FIXME: debugging print
            logger.debug(f"REF={ref}, ALT={alts}, original={entry_allele} -> fixed={fixed_allele}")

    return fixed_entry
