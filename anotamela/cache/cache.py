import json
import logging


logger = logging.getLogger(__name__)


class Cache:
    """
    Base class for caching that provides a common logic to set and retrieve
    group of ids in the form of dictionaries. Particularly, it will deal with
    serializing and deserializing data during dump and load, if needed.

    Classes inheriting from this class should implement the methods:

        - _client_get(ids, namespace)  # => info_dict { id1: val1, id2: val2 ... }
        - _client_set(info_dict, namespace)

    """
    def get(self, ids, namespace, verbose=True):
        """
        Get data for a list of ids from cache. Returns an info dict like:

        {
            id1: raw response of the service for id 1,
            id2: raw response of the service for id 2,
            ...
        }

        During retrieval, the <namespace> will be used as a tablename or id
        prefix according to each particular Cache class.

        Dicts and lists will be json-loaded.
        """
        cached_data = self._client_get(ids, namespace=namespace)
        info_dict = {k: v for k, v in cached_data.items() if v}

        for id_, raw_response in info_dict.items():
            response = self._try_to_decode_and_deserialize(raw_response)
            info_dict[id_] = response

        if verbose and len(ids) > 1:
            msg = 'Found {}/{} ({})'
            logger.info(msg.format(len(info_dict), len(ids),
                                   self.__class__.__name__))

        return info_dict

    def set(self, info_dict, namespace):
        """
        Set the cache for a list of ids. Expects a dict with the form:

        {
            id1: raw response of the service for id 1,
            id2: raw response of the service for id 2,
            ...
        }

        The namespace will be used as a tablename or id prefix according to
        each particular Cache class. Dicts and lists will be json-dumped.
        """
        data_to_cache = {}
        for id_, value in info_dict.items():
            if not value:
                continue
            if isinstance(value, dict) or isinstance(value, list):
                value = json.dumps(value)
            data_to_cache[id_] = value

        if data_to_cache:
            self._client_set(data_to_cache, namespace=namespace)

    @staticmethod
    def _try_to_decode_and_deserialize(data):
        try:
            decoded = data.decode('utf-8')
        except AttributeError:  # Already an UTF-8 str, leave as is
            decoded = data

        try:
            return json.loads(decoded)
        except ValueError:  # Not valid JSON, leave as is
            return decoded

    def _client_get(ids, namespace):
        raise NotImplementedError

    def _client_set(data_to_cache, namespace):
        raise NotImplementedError

