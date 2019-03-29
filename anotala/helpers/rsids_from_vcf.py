def rsids_from_vcf(vcf_path):
    with open(vcf_path) as f:
        geno_lines = [line for line in f if line and not line.startswith('#')]

    ids = [geno_line.split('\t')[2] for geno_line in geno_lines]
    # ^ 3rd position (index=2) in the VCF lines is the ID
    rs_ids = [id_ for id_ in ids if id_.startswith('rs')]

    return rs_ids
