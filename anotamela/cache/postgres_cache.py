from os.path import expanduser
import yaml
import reprlib
import datetime

from sqlalchemy import create_engine, Column, String, DateTime
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker


class PostgresCache():
    def __init__(self, credentials_filepath):
        """
        Accepts a filepath to a YAML file with keys: host, user, pass, port,
        db. It will try to connect to a PostgreSQL database with those
        credentials.
        """
        self._credentials = self._read_credentials()
        self._connect(self._credentials)

        self.models = {}

    def _client_get(self, ids, namespace):
        pass

    def _client_set(self, info_dict, namespace):
        pass

    def get_model(self, tablename):
        """Get a SQLAlchemy model for the given table."""
        if not tablename in self.models:
            self.models[tablename] = self._create_model(tablename)

        return self.models[tablename]

    @staticmethod
    def _connect(credentials):
        url = 'postgresql://{user}:{pass}@{host}:{port}/{db}'.format(
                **credentials)
        self.connection = create_engine(url)
        self.session = sessionmaker().configure(bind=self.connection)

    @staticmethod
    def _read_credentials():
        with open(expanduser('~/.postgres_credentials.yml')) as f:
            return yaml.load(f.read())

    def _create_model(self, tablename):
        """
        Create a SQLAlchemy model for the given tablename. All models created
        have the same columns: id, annotation, synonyms, and last_updated.

        If the table doesn't exist, it will be created in the process.

        Returns the model class, which you can use to query the table.
        """
        def _model_repr(self):
            rep = "{}(id='{}', annotation='{}', synonyms='{}')"
            return rep.format(self.name.capitalize(), self.id,
                            reprlib.repr(self.annotation), self.synonyms)

        class_attributes = {
                '__tablename__': tablename,

                # Table columns:
                'id': Column(String(60), primary_key=True, nullable=False),
                'annotation': Column(JSONB),
                'synonyms': Column(JSONB),
                'last_updated': Column(DateTime(timezone=True), nullable=False,
                                       default=datetime.datetime.utcnow),

                '__repr__': _model_repr,
            }

        Base = declarative_base()
        model_class = type(tablename.capitalize(), (Base, ), class_attributes)

        # This step creates the table if it doesn't exist
        Base.metadata.create_all(self.connection)

        return model_class

