import json

from anotamela.annotators import DbsnpWebAnnotator
from anotamela.helpers import rsids_from_vcf


def annotate_vcf_rsids_with_dbsnp_web(vcf_path, output_json_path=None,
                                      **annotator_options):
    """
    Given a VCF path, annotate its rs IDs with dbSNP Web.

    Returns a dictionary of annotations, or writes the results to a JSON file
    if *output_json_path* is passed.
    """
    rs_ids = rsids_from_vcf(vcf_path)
    annotations = annotate_rsids_with_dbsnp_web(rs_ids, **annotator_options)

    if output_json_path:
        with open(output_json_path, 'w') as f:
            json.dump(annotations, f)
        return output_json_path
    else:
        return annotations



def annotate_rsids_with_dbsnp_web(rs_ids, **annotator_options):
    """
    Annotate a list of rs IDs with dbSNP web.
    Returns a dictionary with the annotations.
    """
    dbsnp_web = DbsnpWebAnnotator(**annotator_options)

    parsed_rs_ids = []
    for id_ in rs_ids:
        if not id_:
            continue

        # Take care of both regular ids like 'rs123', and also of multiple ids
        # like 'rs123;rs234' that sometimes are seen in VCFs:
        parsed_rs_ids += [rs for rs in id_.split(';') if rs.startswith('rs')]

    annotations = dbsnp_web.annotate(parsed_rs_ids)
    return list(annotations.values())
