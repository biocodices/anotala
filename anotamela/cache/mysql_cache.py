from sqlalchemy.dialects.mysql import JSON

from anotamela.cache import SqlCache


class MysqlCache(SqlCache):
    CREDS_FILE = '~/.mysql_credentials.yml'
    JSON_TYPE = JSON

