from collections import defaultdict

import vcf
import pandas as pd
from tqdm import tqdm


def _header_from_vcf(vcf_path):
    """Read the header from a VCF file and return the found field names."""
    with open(vcf_path, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('#CHROM'):
                return line.strip().split('\t')


def vcf_to_dataframe(vcf_path):
    """Generate a pandas.DataFrame of the variants present in a VCF. Removes
    the duplicates."""
    variants = []
    infos = {}

    reader = vcf.Reader(filename=vcf_path)

    for ix, record in enumerate(tqdm(reader)):
        variant = _vcf_record_to_variant_dict(record)
        variant['ix'] = ix  # We use this for latter mapping to INFOs
        infos[ix] = record.INFO
        variants.append(variant)

    df = pd.DataFrame(variants).drop_duplicates()
    df['chrom'] = make_chromosome_series_categorical(df['chrom'])
    df['info'] = df['ix'].map(infos)

    #  for sample, genotype_num_calls in genotypes_nums.items():
    for sample, genotype_num_calls in genotypes_bases.items():
        df[sample] = genotype_num_calls

    variant_cols = 'chrom pos id ref alt qual filter info'.split()
    return df[variant_cols + reader.samples].copy()

def _vcf_record_to_variant_dict(record):
    return {'id': record.ID,
            'chrom': record.CHROM,
            'pos': record.POS,
            'ref': record.REF,
            'qual': record.QUAL,
            'filter': tuple(record.FILTER) or None,
            'alt': tuple(str(allele) for allele in record.ALT)}


