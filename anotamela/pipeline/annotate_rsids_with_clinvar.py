from anotamela.pipeline import annotate_items_with_clinvar


def annotate_rsids_with_clinvar(rs_ids, clinvar_vcf_path=None,
                                grouped_by_rsid=False,
                                **anotator_options):
    """
    Annotate a list of rs IDs with ClinVar (get their ClinVar Variation Report).

    Returns a list of annotations (each one a dictionary).

    Pass a +clinvar_vcf_path+ if you want to specify the path to the ClinVar
    VCF file that will be used to map rs IDs -> variation IDs.

    Pass +grouped_by_rsid+=True to get instead a dictionary with the original rs_ids as
    keys and a list of Clinvar Variation Reports for each of them.
    """
    return annotate_items_with_clinvar(rs_ids, 'rs_id', clinvar_vcf_path,
                                       grouped_by_unique_item=grouped_by_rsid,
                                       **anotator_options)
