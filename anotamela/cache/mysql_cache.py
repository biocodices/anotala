import logging

from sqlalchemy.dialects.mysql import JSON

from anotamela.cache import SqlCache


logger = logging.getLogger(__name__)


class MysqlCache(SqlCache):
    CREDS_FILE = '~/.mysql_credentials.yml'
    JSON_TYPE = JSON

