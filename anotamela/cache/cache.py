import json
import string
import logging


logger = logging.getLogger(__name__)


class Cache:
    """
    Base class for caching that provides a common logic to set and retrieve
    group of ids in the form of dictionaries. Particularly, it will deal with
    serializing and deserializing data during dump and load, if needed.

    Classes inheriting from this class should implement the methods:

        - _client_get(ids, namespace, json_type)
          # => info_dict { id1: val1, id2: val2 ... }
        - _client_set(info_dict, namespace, json_type)

    """
    def get(self, ids, namespace, as_json):
        """
        Get data for a list of ids from cache. Returns an info dict like:

        {
            id1: raw response of the service for id 1,
            id2: raw response of the service for id 2,
            ...
        }

        During retrieval, the <namespace> will be used as a tablename or id
        prefix according to each particular Cache class.

        If as_json = True, values will be json-loaded.
        """
        cached_data = self._client_get(ids, namespace, as_json)

        logger.info('{} found {}/{}'.format(self.__class__.__name__,
                                            len(cached_data), len(ids)))
        return cached_data

    def set(self, info_dict, namespace, as_json):
        """
        Set the cache for a list of ids. Expects a dict with the form:

        {
            id1: raw response of the service for id 1,
            id2: raw response of the service for id 2,
            ...
        }

        The namespace will be used as a tablename or id prefix according to
        each particular Cache class.

        If as_json = True, the values will be json-dumped.
        """
        data_for_cache = {}

        for id_, annotation in info_dict.items():
            if annotation:

                if not as_json:
                    annotation = self._only_printable(annotation)

            data_for_cache[id_] = annotation

        if data_for_cache:
            self._client_set(data_for_cache, namespace, as_json)

    @staticmethod
    def _only_printable(s):
        return ''.join(char for char in s if char in string.printable)

    @staticmethod
    def _jsondump_dict_values(dictionary):
        return {key: json.dumps(value) for key, value in dictionary.items()}

    @staticmethod
    def _jsonload_dict_values(dictionary):
        return {key: json.loads(value) for key, value in dictionary.items()}

    def _client_get(self):
        raise NotImplementedError

    def _client_set(self):
        raise NotImplementedError

