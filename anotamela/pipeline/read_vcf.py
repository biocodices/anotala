def read_variants_from_vcf(vcf_path):
    """
    Read the VCF at *vcf_path* and generate pandas DataFrames for:

        - the variants with rs ID
        - the other variants (i.e. without rs ID)

    Returns a dictionary with both dataframes under the keys 'rs_variants' and
    'other_variants'.
    """
    logger.info('Read "{}"'.format(vcf_path))
    self.variants = vcf_to_dataframe(vcf_path)
    have_single_rs = self.variants['id'].str.match(r'^rs\d+$')
    self.rs_variants = self.variants[have_single_rs].reset_index(drop=True)
    self.other_variants = self.variants[~have_single_rs].reset_index(drop=True)
    logger.info('{} variants with single rs'.format(len(self.rs_variants)))
    logger.info('{} other variants'.format(len(self.other_variants)))


