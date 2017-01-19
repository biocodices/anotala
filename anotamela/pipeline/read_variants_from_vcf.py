from vcf_to_dataframe import vcf_to_dataframe


def read_variants_from_vcf(vcf_path):
    """
    Read the VCF at *vcf_path* and generate pandas DataFrames for:

        - the variants with rs ID
        - the other variants (i.e. without rs ID)

    Returns a dictionary with both dataframes under the keys 'rs_variants' and
    'other_variants'.
    """
    variants = vcf_to_dataframe(vcf_path)
    have_single_rs = variants['id'].str.match(r'^rs\d+$')

    return {
        'rs_variants': variants[have_single_rs].reset_index(drop=True),
        'other_variants': variants[~have_single_rs].reset_index(drop=True),
    }

