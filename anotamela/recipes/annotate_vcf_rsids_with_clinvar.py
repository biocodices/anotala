import json

from anotamela.helpers import rsids_from_vcf
from anotamela.pipeline import annotate_rsids_with_clinvar


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
