def generate_position_tags(variant, assembly):
    """
    Give a variant with a 'dbsnp_web' key, return a position tag like "1:10000"
    according to the chosen assembly.
    """
    chrom_key = f'{assembly}_chrom'
    pos_key = f'{assembly}_start'
    chromosome = variant['dbsnp_web'].get(chrom_key) or ''
    position = variant['dbsnp_web'].get(pos_key) or ''
    return f'{chromosome}:{position}'
