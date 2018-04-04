import json
from collections import defaultdict

from anotamela.annotators import (
    ClinvarVariationAnnotator,
    ClinvarRsVCFAnnotator,
)
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


def _get_clinvar_variation_ids_from_rs_ids(rs_ids, clinvar_vcf_path):
    clinvar_vcf = ClinvarRsVCFAnnotator(path_to_annotations_file=clinvar_vcf_path)
    annotations = clinvar_vcf.annotate(rs_ids)
    variation_ids = set()
    for rs, vcf_entries in annotations.items():
        for vcf_entry in vcf_entries:
            variation_ids.add(vcf_entry['variation_id'])
    return variation_ids


def annotate_rsids_with_clinvar(rs_ids, cache, clinvar_vcf_path=None,
                                proxies={}, use_web=True, use_cache=True,
                                grouped_by_rsid=False):
    """
    Annotate a list of rs IDs with ClinVar (get their ClinVar Variation Report).

    Returns a list of annotations (each one a dictionary).

    Pass a +clinvar_vcf_path+ if you want to specify the pat to the ClinVar
    VCF file that will be used to map rs IDs -> variation IDs.

    Pass +grouped_by_rsid+=True to get instead a dictionary with the original rs_ids as
    keys and a list of Clinvar Variation Reports for each of them.
    """
    annotator_options = dict(cache=cache, proxies=proxies)
    annotate_options = dict(use_web=use_web, use_cache=use_cache)

    variation_ids = \
        _get_clinvar_variation_ids_from_rs_ids(rs_ids,
                                               clinvar_vcf_path=clinvar_vcf_path)

    clinvar_var = ClinvarVariationAnnotator(**annotator_options)
    clinvar_variations = clinvar_var.annotate(variation_ids, **annotate_options)
    clinvar_variations = list(clinvar_variations.values())

    if grouped_by_rsid:
        clinvar_variations_by_rsid = defaultdict(list)

        for variation in clinvar_variations:
            dbsnp_id = variation.get('dbsnp_id')
            dbsnp_ids = variation.get('dbsnp_ids')

            if dbsnp_id:
                clinvar_variations_by_rsid[dbsnp_id].append(variation)
            elif dbsnp_ids:
                for id_ in dbsnp_ids:
                    clinvar_variations_by_rsid[id_].append(variation)

        # Make sure all queried rs_ids are present in the final dictionary:
        for rs_id in rs_ids:
            if rs_id not in clinvar_variations_by_rsid:
                clinvar_variations_by_rsid[rs_id] = []

        # Remove clinvar reports about other variants:
        for key in list(clinvar_variations_by_rsid.keys()):
            if key not in rs_ids:
                del(clinvar_variations_by_rsid[key])

        clinvar_variations = clinvar_variations_by_rsid

    return clinvar_variations
