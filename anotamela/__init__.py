from os import environ
import logging

import coloredlogs

from .annotators import DbSNPAnnotator


coloredlogs.DEFAULT_LOG_FORMAT = \
        '[@%(hostname)s %(asctime)s] %(levelname)s: %(message)s'
coloredlogs.install(level=environ.get('LOGLEVEL') or 'INFO')

logging.getLogger('requests').setLevel(logging.WARNING)
logging.getLogger('urllib3').setLevel(logging.WARNING)

