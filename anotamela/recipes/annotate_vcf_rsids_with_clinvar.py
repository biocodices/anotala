import json

from anotamela.annotators import ClinvarRsAnnotator, ClinvarVariationAnnotator
from anotamela.helpers import rsids_from_vcf


def annotate_vcf_rsids_with_clinvar(vcf_path, output_json_path=None,
                                    **annotator_options):
    """
    Given a VCF path, annotate its rs IDs with Clinvar.

    Returns a dictionary of annotations, or writes the results to a JSON file
    if *output_json_path* is passed.
    """
    rs_ids = rsids_from_vcf(vcf_path)
    annotations = annotate_rsids_with_clinvar(rs_ids, **annotator_options)

    if output_json_path:
        with open(output_json_path, 'w') as f:
            json.dump(annotations, f)
        return output_json_path
    else:
        return annotations


def annotate_rsids_with_clinvar(rs_ids, **annotator_options):
    """
    Annotate a list of rs IDs with ClinVar.
    Returns a dictionary with the annotations.
    """
    clinvar_rs = ClinvarRsAnnotator(**annotator_options)
    annotations_rs = clinvar_rs.annotate(rs_ids)

    variant_ids = [entry.get('variant_id')
                   for annotation in annotations_rs.values()
                   for entry in annotation]
    variant_ids = [id_ for id_ in variant_ids if id_]

    clinvar_var = ClinvarVariationAnnotator(**annotator_options)
    annotations = clinvar_var.annotate(variant_ids)

    return list(annotations.values())
