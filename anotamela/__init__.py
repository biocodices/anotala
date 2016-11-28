from os import environ
import coloredlogs


coloredlogs.install(level=environ.get('LOGLEVEL') or 'INFO')
coloredlogs.DEFAULT_LOG_FORMAT = '[@%(hostname)s %(asctime)s] %(name)s : %(message)s'

