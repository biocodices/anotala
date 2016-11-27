from collections import defaultdict
import logging

import coloredlogs

from anotamela.helpers import make_chromosome_series_categorical
from anotamela.vcf_to_dataframe import vcf_to_dataframe


logger = logging.getLogger(__name__)
coloredlogs.DEFAULT_LOG_FORMAT = '[@%(hostname)s %(asctime)s] %(message)s'
coloredlogs.install(level='INFO')

def pipeline(vcf_path):
    """
    Annotate the given VCF file. Return the filepath to the annotated (new)
    VCF.
    """
    logger.info('Read "{}"'.format(vcf_path))
    variants = vcf_to_dataframe(vcf_path)
    return variants

