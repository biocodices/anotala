from anotamela.cache import AVAILABLE_CACHES


AVAILABLE_CACHES['mock_cache'] = AVAILABLE_CACHES['_dict']
PROXIES = {'http': 'socks5://beleriand.local:9150'}

