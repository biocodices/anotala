import re

import vcf
import gzip
import pandas as pd
from tqdm import tqdm

from anotamela.helpers import make_chromosome_series_categorical


GENO_REGEX = re.compile(r'([\d|\.](?:[/|\|][\d|\.])?)')

def vcf_to_dataframe(vcf_path, keep_format_data=False, keep_samples=[]):
    """
    Generate a pandas.DataFrame of the variants present in a VCF. Removes
    the duplicated rows. WARNING: It will hog RAM and take a while for big
    VCFs.

    To avoid big time use of RAM and cycles, the default behavior is not to
    read any of the genotypes. In case you want to keep the genotypes, set
    keep_samples with a sample ID or list of sample IDs. In case you want
    all the samples, set keep_samples='all' (WARNING: this will take loads of
    RAM and time for big VCFs).

    If keep_format_data is set, it will keep the metadata for each
    genotype call, e.g. AD, DP, GQ, PL, etc. If not, just the genotypes (GT).
    """
    header = _header_from_vcf(vcf_path, gzipped=vcf_path.endswith('.gz'))

    if keep_samples:
        if isinstance(keep_samples, str):
            keep_samples = [keep_samples]

        for sample in keep_samples:
            if not sample in header:
                raise ValueError('"{}" not found in this VCF'.format(sample))

    usecols = header[:9] + keep_samples

    df = pd.read_table(vcf_path, comment='#', names=header, sep='\s+',
                       low_memory=False, usecols=usecols)
    df.rename(columns={'#CHROM': 'CHROM'}, inplace=True)
    df.rename(columns={col: col.lower() for col in df.columns[:9]}, inplace=True)

    if keep_samples and not keep_format_data:
        _extract_genos_and_make_them_categorical(df)

    for col in 'ref alt filter'.split():
        df[col] = df[col].astype('category')

    df['chrom'] = make_chromosome_series_categorical(df['chrom'])

    return df

def _extract_genos_and_make_them_categorical(df):
    """Extract the genotypes from the format fields of the VCF and make them
    categorical. It changes the input DataFrame!"""
    # Genotype columns range from the 10th to the last one
    df.iloc[:, 9:] = df.iloc[:, 9:].applymap(_extract_genotype)

    # Cast genotypes as category since it's MUCH more memory efficient
    df.iloc[:, 9:] = df.iloc[:, 9:].apply(
            lambda series: series.astype('category'))

def _header_from_vcf(vcf_path, gzipped=False):
    """Read the header from a VCF file and return the found field names. It
    will gunzip on the fly if the path ends with '.gz' or if gzipped is set."""
    fn_open = gzip.open if gzipped else open

    with fn_open(vcf_path, 'rb') as vcf_file:
        for line in vcf_file:
            if isinstance(line, bytes):
                line = line.decode('utf-8')
            if line.startswith('#CHROM'):
                return line.strip().split('\t')

def _extract_genotype(geno_field):
    """Extract the genotype from a format field."""
    # Assume the genotype is the first format field, raise if it's not.
    geno = geno_field.split(':')[0]
    if not GENO_REGEX.search(geno):
        raise ValueError('"{}" does not look like a genotype'.format(geno))
    return geno

