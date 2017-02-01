from sqlalchemy.dialects.mysql import JSON, MEDIUMTEXT

from anotamela.cache import SqlCache


class MysqlCache(SqlCache):
    CREDS_FILE = '~/.mysql_credentials.yml'
    JSON_TYPE = JSON
    TEXT_TYPE = MEDIUMTEXT
    URL = '{driver}://{user}:{pass}@{host}:{port}/{db}?charset=utf8'

