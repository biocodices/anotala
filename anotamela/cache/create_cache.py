from anotamela.cache import PostgresCache, RedisCache


AVAILABLE_CACHES = {
        'redis': RedisCache,
        'postgres': PostgresCache
    }


def create_cache(cache_name, **cache_kwargs):
    try:
        cache_class = AVAILABLE_CACHES[cache_name]
    except KeyError as error:
        known_caches = ', '.join(AVAILABLE_CACHES.keys())
        msg = 'Unknown cache "{}". I only know: {}'.format(cache_name,
                                                           known_caches)
        raise ValueError(msg).with_traceback(error.__traceback__)
    else:
        return cache_class(**cache_kwargs)

