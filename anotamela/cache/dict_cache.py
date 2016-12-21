import json
from collections import defaultdict

from anotamela.cache import Cache


class DictCache(Cache):
    """
    This class is used for testing and housekeeping scripts. It shouldn't be
    used for the actual variant annotation pipelines, since it's not thread
    safe and it won't escalate properly.
    """
    def __init__(self):
        self.name = 'DictCache'
        self.storage = defaultdict(dict)

    def _client_set(self, info_dict, namespace, save_as_json):
        table = self.storage[namespace]
        if save_as_json:
            info_dict = {k: json.dumps(v) for k, v in info_dict.items()}
        table.update(info_dict)

    def _client_get(self, keys, namespace, load_as_json):
        table = self.storage[namespace]
        info_dict = {k: table[k] for k in keys if k in table}
        if load_as_json:
            info_dict = {k: json.loads(v) for k, v in info_dict.items()}
        return info_dict

