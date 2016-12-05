from os.path import expanduser
import yaml
import logging
from datetime import datetime

from sqlalchemy import create_engine, Table, Column, String, DateTime, MetaData
from sqlalchemy.dialects.postgresql import JSONB

from anotamela.cache import Cache


logger = logging.getLogger(__name__)


class PostgresCache(Cache):
    CREDS_FILE = '~/.postgres_credentials.yml'
    SAVE_AS_JSON = False  # Don't assume JSON types unless it's explicitely set

    def __init__(self, credentials_filepath=CREDS_FILE):
        """
        Accepts a filepath to a YAML file with keys: host, user, pass, port,
        db. It will try to connect to a PostgreSQL database with those
        credentials.
        """
        self._credentials = self._read_credentials(credentials_filepath)
        self._connect(self._credentials)
        logger.info('Connected to PostgreSQL ({!r})'.format(self.engine.url))
        self.tables = {}

    def _client_get(self, ids, namespace):
        table = self._get_table(namespace)
        selection = table.select().where(table.c.id.in_(ids))
        result = self.connection.execute(selection)
        return {row['id']: row['annotation'] for row in result}

    def _client_set(self, info_dict, namespace):
        table = self._get_table(namespace)
        ids_to_remove = info_dict.keys()
        new_annotations = [{'id': id_, 'annotation': ann}
                           for id_, ann in info_dict.items()]

        self._client_del(ids_to_remove, namespace)
        self.connection.execute(table.insert(), new_annotations)

    def _client_del(self, ids, namespace):
        table = self._get_table(namespace)
        deletion = table.delete().where(table.c.id.in_(ids))
        self.connection.execute(deletion)

    def _connect(self, credentials):
        url_template = 'postgresql://{user}:{pass}@{host}:{port}/{db}'
        url = url_template.format(**credentials)
        self.engine = create_engine(url)
        self.connection = self.engine.connect()

    def _get_table(self, tablename):
        """Get a SQLAlchemy Table for the given tablename."""
        if tablename not in self.tables:
            self.tables[tablename] = self._create_table(tablename)

        return self.tables[tablename]

    def _create_table(self, tablename):
        """
        Create a SQLAlchemy Table for the given tablename. All tables created
        have the same columns: id, annotation, synonyms, and last_updated.

        If the table doesn't exist already in the database, it will be created.

        Returns a sqlalchemy.Table instance.
        """
        metadata = MetaData()
        annotation_column_type = JSONB if self.SAVE_AS_JSON else String
        table = Table(tablename, metadata,
                      Column('id', String(60), primary_key=True),
                      Column('annotation', annotation_column_type),
                      Column('synonyms', JSONB),  # Just a list of IDs
                      Column('last_updated', DateTime(timezone=True),
                             default=datetime.now, nullable=False))
        metadata.create_all(self.engine)
        return table

    @staticmethod
    def _read_credentials(filepath):
        try:
            with open(expanduser(filepath)) as f:
                return yaml.load(f.read())
        except FileNotFoundError as error:
            msg = "Couldn't find a YAML with PostgreSQL credentials in '{}'"
            msg += '. It should include: host, user, pass, port, db.'
            msg = msg.format(filepath)
            raise FileNotFoundError(msg).with_traceback(error.__traceback__)

