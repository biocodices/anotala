import json
import logging


logger = logging.getLogger(__name__)


class Cache:
    """
    Base class for caching that provides a common logic to set and retrieve
    group of keys in the form of info dictionaries. Particularly, it will deal
    with serializing and deserializing data during dump and load.

    Classes inheriting from this class should implement the methods:

        - _client_get(keys)
        - _client_set(info_dict)

    """
    def _client_get(keys):
        raise NotImplementedError

    def _client_set(info_dict):
        raise NotImplementedError

    def get(self, keys, verbose=True):
        """
        Get a list of keys data from cache. Returns an info dict with the form:

        {
            key1: raw response of the service for key 1,
            key2: raw response of the service for key 2,
            ...
        }

        It will json-load the dicts and lists.
        """
        info_dict = {k: v for k, v in zip(keys, self._client_get(keys)) if v}

        for key, raw_response in info_dict.items():
            response = self._decode_and_try_to_deserialize(raw_response)
            info_dict[key] = response

        if verbose and len(keys) > 1:
            msg = '📂 {} found {}/{} in cache'
            logger.info(msg.format(self.name, len(info_dict), len(keys)))

        return info_dict

    def set(self, info_dict):
        """
        Set the cache for a list of keys. Expects a dict with the form:

        {
            key1: raw response of the service for key 1,
            key2: raw response of the service for key 2,
            ...
        }

        It will json-dump the dicts and lists.
        """
        data_to_cache = {}
        for key, value in info_dict.items():
            if not value:
                continue
            if isinstance(value, dict) or isinstance(value, list):
                value = json.dumps(value)
            data_to_cache[key] = value

        self._client_set(data_to_cache)

    @staticmethod
    def _decode_and_try_to_deserialize(data):
        try:
            decoded = data.decode('utf-8')
        except AttributeError:  # Already an UTF-8 str, leave as is
            decoded = data

        try:
            return json.loads(decoded)
        except ValueError:  # Not valid JSON, leave as is
            return decoded

