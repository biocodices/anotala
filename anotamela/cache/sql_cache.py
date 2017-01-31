from os.path import expanduser, abspath
import yaml
import logging

from sqlalchemy import Table, Column, String, DateTime, MetaData, func

from anotamela.cache import Cache


logger = logging.getLogger(__name__)


class SqlCache(Cache):
    """
    Base class for PostgresCache and MysqlCache. To implement a SqlCache class,
    define:

        - CREDS_FILE filepath with credentials for the database, in yaml.
          For instance:

            host: localhost
            port: 5432
            db: variants_cache
            user: johndoe
            pass: jellybeans

            driver: (mysql+pymysql|postgresql)
            ^ See : http://docs.sqlalchemy.org/en/latest/core/engines.html

        - JSON_TYPE for SQLAlchemy to use, e.g.

          - sqlalchemy.dialects.postgresql.JSONB for PostgreSQL
          - sqlalchemy.dialects.mysql.JSON for MySQL

    Provided a database has already been created and its name and credentials
    are stored in the *credentials_filepath*, this class will connect with the
    db and deal with table creation and cache get/set.

    Tables in the database will be created using the passed namespaces as names
    and they will always have the same fields: 'id', 'annotation', 'synonyms',
    and a 'last_updated' timestamp.
    """

    CREDS_FILE = None  # Overwrite this value in the subclass
    JSON_TYPE = None  # Overwrite this value in the subclass
    DIALECT_DRIVER = None  # Overwrite this value in the subclass

    def __init__(self, credentials_filepath=None):
        """
        Accepts a filepath to a YAML file with keys: host, user, pass, port,
        db. It will try to connect to a PostgreSQL database with those
        credentials.
        """
        self.credentials_filepath = \
            abspath(expanduser(credentials_filepath or self.CREDS_FILE))
        self.credentials = self._read_credentials(self.credentials_filepath)
        self.metadata = self._connect_to_db(self.credentials)
        self.engine = self.metadata.bind  # Shortcut for statements execution

        logger.info('Connected to {} ({!r})'
                    .format(self.__class__.__name__, self.engine.url))

    def __repr__(self):
        return "{}('{}')".format(self.__class__.__name__,
                                 self.credentials_filepath)

    def _client_get(self, ids, namespace, as_json):
        """
        Gets data for a list of IDs in the given namespace.

        If as_json is True, the table will be created with JSON type the
        first time a namespace is used, and in later uses the data will be
        deserialized as JSON.

        Returns a dictionary with the IDs as keys and their cache as values.
        """
        # The argument as_json is only relevant the first time an
        # annotator is used, in self._get_table(), to create the table with a
        # JSON field type. After this, the SQLAlchemy automatically deals with
        # [de]serialization of JSON because the fields are set as JSON type.
        table = self._get_table(namespace, as_json=as_json)

        select_query = table.select().where(table.c.id.in_(ids))
        query_result = self.engine.execute(select_query)
        return {row['id']: row['annotation'] for row in query_result}

    def _client_set(self, info_dict, namespace, as_json):
        """
        Accepts and *info_dict* with IDs as keys and annotations as values.
        It will write the annotations to a table named after the *namespace*
        in the database. If some IDs already exist, those rows will be
        overwritten.
        """
        # The argument as_json is only relevant the first time an
        # annotator is used, in self._get_table(), to create the table with a
        # JSON field type in case the table named *namespace* doesn't exist.
        table = self._get_table(namespace, as_json=as_json)

        ids_to_remove = info_dict.keys()
        self._client_del(ids_to_remove, namespace)
        new_annotations = [{'id': id_, 'annotation': ann}
                           for id_, ann in info_dict.items()]
        self.engine.execute(table.insert(), new_annotations)

    def _client_del(self, ids, namespace):
        """
        Delete the passed *ids* from the table *namespace* in the database.
        """
        table = self._get_table(namespace)
        delete_query = table.delete().where(table.c.id.in_(ids))
        self.engine.execute(delete_query)

    @staticmethod
    def _connect_to_db(credentials):
        """
        Connect to the database and read the available tables.
        Returns a SQLAlchemy MetaData instance.
        """
        url = '{driver}://{user}:{pass}@{host}:{port}/{db}'
        metadata = MetaData(bind=url.format(**credentials))
        metadata.reflect()  # Gets a list of all Tables in the database
        return metadata

    def _get_table(self, tablename, as_json=False):
        """
        Get a SQLAlchemy Table for the given tablename. If it doesn't exist
        in the database, it will be created.
        """
        if tablename not in self.metadata.tables:
            self._create_table(tablename, as_json)

        return self.metadata.tables[tablename]

    def _create_table(self, tablename, as_json):
        """
        Create a SQLAlchemy Table for the given tablename. All tables created
        have the same columns: id, annotation, synonyms, and last_updated.

        If the table doesn't exist already in the database, it will be created.

        Returns a sqlalchemy.Table instance.
        """
        print('Create {} called'.format(tablename))
        annotation_column_type = self.JSON_TYPE if as_json else String

        Table(tablename, self.metadata,
              Column('id', String(60), primary_key=True),
              Column('annotation', annotation_column_type),
              Column('synonyms', self.JSON_TYPE),  # Just a list of IDs
              Column('last_updated', DateTime(timezone=True),
                     server_default=func.now(), nullable=False))

        self.metadata.create_all(checkfirst=True)

        return self.metadata.tables[tablename]

    @staticmethod
    def _read_credentials(filepath):
        """
        Read credentials to access the database from the given filepath. It
        should point to a YAML file with: host, user, pass, port, database,
        driver. Example:

            host: localhost
            port: 3306
            db: variants_cache
            user: johndoe
            pass: jellybeans
            driver: mysql

        """
        try:
            with open(expanduser(filepath)) as f:
                return yaml.load(f.read())
        except IOError as error:
            msg = "Couldn't read a YAML with PostgreSQL credentials in '{}'"
            msg += '. It should include: host, user, pass, port, db.'
            msg = msg.format(filepath)
            raise IOError(msg).with_traceback(error.__traceback__)

