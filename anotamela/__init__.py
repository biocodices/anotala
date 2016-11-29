from os import environ
import coloredlogs

from .annotators import DbSNPAnnotator


coloredlogs.DEFAULT_LOG_FORMAT = \
        '[@%(hostname)s %(asctime)s] %(levelname)s: %(message)s'
coloredlogs.install(level=environ.get('LOGLEVEL') or 'INFO')

