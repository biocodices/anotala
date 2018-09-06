from anotamela.pipeline import annotate_items_with_clinvar


def annotate_position_tags_with_clinvar(position_tags, clinvar_vcf_path=None,
                                        grouped_by_position=False,
                                        **anotator_options):
    """
    Annotate a list of position_tags like "1:10000" with ClinVar (get their
    ClinVar Variation Reports). Returns a list of annotations (each one a
    dictionary).

    Pass a +clinvar_vcf_path+ if you want to specify the pat to the ClinVar
    VCF file that will be used to map rs IDs -> variation IDs.

    Pass grouped_by_position=True to get instead a dictionary with the original
    positions as keys and a list of Clinvar Variation Reports for each of them.
    """
    return annotate_items_with_clinvar(
        position_tags, 'chrom_pos', clinvar_vcf_path,
        grouped_by_unique_item=grouped_by_position, **anotator_options
    )



