import logging

from sqlalchemy.dialects.postgresql import JSONB

from anotamela.cache import SqlCache


logger = logging.getLogger(__name__)


class PostgresCache(SqlCache):
    CREDS_FILE = '~/.postgres_credentials.yml'
    JSON_TYPE = JSONB

