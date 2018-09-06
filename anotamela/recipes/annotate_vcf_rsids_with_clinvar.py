import json
from collections import defaultdict

from more_itertools import collapse, unique_everseen

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


def annotate_rsids_with_clinvar(rs_ids, clinvar_vcf_path=None,
                                grouped_by_rsid=False,
                                **anotator_options):
    """
    Annotate a list of rs IDs with ClinVar (get their ClinVar Variation Report).

    Returns a list of annotations (each one a dictionary).

    Pass a +clinvar_vcf_path+ if you want to specify the pat to the ClinVar
    VCF file that will be used to map rs IDs -> variation IDs.

    Pass +grouped_by_rsid+=True to get instead a dictionary with the original rs_ids as
    keys and a list of Clinvar Variation Reports for each of them.
    """
    return annotate_items_with_clinvar(rs_ids, 'rs_id', clinvar_vcf_path,
                                       grouped_by_unique_item=grouped_by_rsid,
                                       **anotator_options)


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


def annotate_items_with_clinvar(items, items_type, clinvar_vcf_path,
                                grouped_by_unique_item=False, cache='mysql',
                                proxies={}, use_web=True, use_cache=True):
    """
    Given a list of *items* like rs IDs and the *items_type*, "rs_id",
    return the Clinvar Variation Reports associated with them.
    """
    items = list(unique_everseen(items))

    # Associate each item (an RS ID, a genomic location, etc.) to a list of
    # one or many ClinVar Variation IDs
    variation_ids_per_item = _get_clinvar_variation_ids_from(
        items, items_type, clinvar_vcf_path=clinvar_vcf_path
    )
    variation_ids = list(collapse(list(variation_ids_per_item.values())))

    # Use the ClinVar Variation IDs to get Variation Reports
    annotator_options = dict(cache=cache, proxies=proxies)
    annotate_options = dict(use_web=use_web, use_cache=use_cache)
    annotator = ClinvarVariationAnnotator(**annotator_options)
    clinvar_variations = annotator.annotate(variation_ids, **annotate_options)

    # Associate the original items with the obtained Variation Reports
    if grouped_by_unique_item:
        result = defaultdict(list)

        for item, item_variation_ids in variation_ids_per_item.items():
            result[item] = [clinvar_variations[variation_id]
                            for variation_id in item_variation_ids
                            if variation_id in clinvar_variations]
    else:
        result = list(clinvar_variations.values())

    return result


def _choose_annotator_class(items_type):
    if items_type == 'chrom_pos':
        annotator_class = ClinvarPosVCFAnnotator
    elif items_type == 'rs_id':
        annotator_class = ClinvarRsVCFAnnotator
    else:
        raise f"I don't know how to filter the ClinVar VCF by {items_type}"
    return annotator_class


def _get_clinvar_variation_ids_from(items, items_type, clinvar_vcf_path):
    annotator_class = _choose_annotator_class(items_type)
    clinvar_vcf_annotator = annotator_class(path_to_annotations_file=clinvar_vcf_path)
    annotations = clinvar_vcf_annotator.annotate(items)
    variation_ids = defaultdict(list)
    for item, clinvar_vcf_entries in annotations.items():
        for vcf_entry in clinvar_vcf_entries:
            # Make the IDs strings so that they match the keys of their
            # annotation downstream:
            variation_id = str(vcf_entry['variation_id'])
            variation_ids[item].append(variation_id)
    return variation_ids
