import logging
import coloredlogs

from vcf_to_dataframe import vcf_to_dataframe
from anotamela.annotators import (
        DbsnpWebAnnotator,
        DbsnpEntrezAnnotator,
        HgvsAnnotator,
        ClinvarRsAnnotator
    )


logger = logging.getLogger(__name__)
coloredlogs.DEFAULT_LOG_FORMAT = '[@%(hostname)s %(asctime)s] %(message)s'
coloredlogs.install(level='INFO')


def pipeline(vcf_path, keep_samples=None, use_cache=True, use_web=True):
    """Annotate the given VCF file. Return the filepath to the annotated (new)
    VCF."""
    logger.info('Read "{}"'.format(vcf_path))
    variants = vcf_to_dataframe(vcf_path, keep_samples=keep_samples)

    single_rs = variants['id'].str.match(r'^rs\d+$')
    rs_variants = variants[single_rs].reset_index(drop=True)
    other_variants = variants[~single_rs].reset_index(drop=True)
    logger.info('%s variants with a single rs ID' % len(rs_variants))
    logger.info("%s other variants won't be annotated" % len(other_variants))

    annotators = [DbsnpWebAnnotator,
                  DbsnpEntrezAnnotator,
                  HgvsAnnotator,
                  ClinvarRsAnnotator]
    for annotator_class in annotators:
        annotator = annotator_class(cache='postgres')
        annotations = annotator.annotate(rs_variants['id'],
                                         use_cache=use_cache, use_web=use_web)
        rs_variants[annotator_class.SOURCE_NAME] = \
                rs_variants['id'].map(annotations)

    return rs_variants

