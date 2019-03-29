from os import environ

import logging

logger = logging.getLogger('anotala')

from Bio import Entrez
import coloredlogs

from .annotators import *
from .pipeline import AnnotationPipeline


coloredlogs.DEFAULT_LOG_FORMAT = \
        '[@%(hostname)s %(asctime)s] %(levelname)s: %(message)s'
coloredlogs.install(level=(environ.get('LOGLEVEL') or 'INFO'), logger=logger)

logging.getLogger('requests').setLevel(logging.WARNING)
logging.getLogger('urllib3').setLevel(logging.WARNING)

Entrez.tool = 'anotala'
