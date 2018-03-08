import json
from itertools import groupby
from operator import itemgetter

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


def annotate_rsids_with_clinvar(rs_ids, cache, proxies={}, use_web=True,
                                use_cache=True, grouped_by_rsid=False):
    """
    Annotate a list of rs IDs with ClinVar (get their ClinVar Variation Report).

    Returns a list of annotations (each one a dictionary).

    Pass +grouped_by_rsid+=True to get instead a dictionary with the original rs_ids as
    keys and a list of Clinvar Variation Reports for each of them.
    """
    annotator_options = dict(cache=cache, proxies=proxies)
    annotate_options = dict(use_web=use_web, use_cache=use_cache)


    clinvar_rs = ClinvarRsAnnotator(**annotator_options)
    annotations_rs = clinvar_rs.annotate(rs_ids, **annotate_options)

    variant_ids = [entry.get('variant_id')
                   for annotation in annotations_rs.values()
                   for entry in annotation]
    variant_ids = [id_ for id_ in variant_ids if id_]

    clinvar_var = ClinvarVariationAnnotator(**annotator_options)
    clinvar_reports = clinvar_var.annotate(variant_ids, **annotate_options)
    clinvar_reports = list(clinvar_reports.values())

    if grouped_by_rsid:
        clinvar_reports = {
            k: list(group) for k, group in
            groupby(clinvar_reports, key=itemgetter('dbsnp_id'))
        }
        # Make sure all queried rs_ids are present in the final dictionary:
        for rs_id in rs_ids:
            if rs_id not in clinvar_reports:
                clinvar_reports[rs_id] = []

    return clinvar_reports
