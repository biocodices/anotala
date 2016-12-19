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

    def __init__(self, credentials_filepath=CREDS_FILE):
        """
        Accepts a filepath to a YAML file with keys: host, user, pass, port,
        db. It will try to connect to a PostgreSQL database with those
        credentials.
        """
        self._creds_filepath = credentials_filepath
        self._credentials = self._read_credentials(credentials_filepath)
        self._connect(self._credentials)
        logger.info('Connected to PostgreSQL ({!r})'.format(self.engine.url))
        self.tables = {}

    def __repr__(self):
        rep = "{name}('{cred}')"
        return rep.format(cred=self._creds_filepath,
                          name=self.__class__.__name__)

    def _client_get(self, ids, namespace, load_as_json):
        # The argument load_as_json is only relevant the first time an
        # annotator is used, in self._get_table(), to create the table with a
        # JSONB field type. After this, the psql client autoatically deals with
        # [de]serialization of JSON, so the setting and retrieval of data is
        # the same wether the field is JSON or not.
        table = self._get_table(namespace, json_type=load_as_json)
        select_query = table.select().where(table.c.id.in_(ids))
        query_result = self.connection.execute(select_query)
        return {row['id']: row['annotation'] for row in query_result}

    def _client_set(self, info_dict, namespace, save_as_json):
        # The argument load_as_json is only relevant the first time an
        # annotator is used, in self._get_table(), to create the table with a
        # JSONB field type. After this, the psql client autoatically deals with
        # [de]serialization of JSON, so the setting and retrieval of data is
        # the same wether the field is JSON or not.
        table = self._get_table(namespace, json_type=save_as_json)
        ids_to_remove = info_dict.keys()
        self._client_del(ids_to_remove, namespace)
        new_annotations = [{'id': id_, 'annotation': ann}
                           for id_, ann in info_dict.items()]
        self.connection.execute(table.insert(), new_annotations)

    def _client_del(self, ids, namespace):
        table = self._get_table(namespace)
        delete_query = table.delete().where(table.c.id.in_(ids))
        self.connection.execute(delete_query)

    def _connect(self, credentials):
        url_template = 'postgresql://{user}:{pass}@{host}:{port}/{db}'
        url = url_template.format(**credentials)
        self.engine = create_engine(url)
        self.connection = self.engine.connect()

    def _get_table(self, tablename, json_type=False):
        """Get a SQLAlchemy Table for the given tablename."""
        if tablename not in self.tables:
            self.tables[tablename] = self._create_table(tablename, json_type)

        return self.tables[tablename]

    def _create_table(self, tablename, json_type):
        """
        Create a SQLAlchemy Table for the given tablename. All tables created
        have the same columns: id, annotation, synonyms, and last_updated.

        If the table doesn't exist already in the database, it will be created.

        Returns a sqlalchemy.Table instance.
        """
        metadata = MetaData()
        annotation_column_type = JSONB if json_type else String
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

