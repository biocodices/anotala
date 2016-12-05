## Annotation

## Usage
```python
from anotamela import (DbsnpMyvariantAnnotator,
                       ClinvarRsAnnotator,
                       OmimVariantAnnotator)

dbsnp_annotator = DbsnpMyvariantAnnotator(cache='redis')  # or cache='postgres'
dbsnp_annotator.annotate(['rs268', 'rs123'])
# => {'rs268': { ... }, 'rs123': { ... }}

omim_annotator = OmimVariantAnnotator(cache='postgres')
dbsnp_entrez.annotate(['605557.0001', '605557.0003'])
# => {'605557.0001': { ... }, '605557.003': { ... }}

# Just fetch from cache, don't use the Internets:
dbsnp_annotator.annotate(rs_list, use_web=False)

# Fetch everything from the Net, don't use the cache.
# This will also update the existing cache entries:
dbsnp_annotator.annotate(rs_list, use_cache=False)
```

There are plenty of annotators: `DbsnpWebAnnotator`, `DbsnpEntrezAnnotator`,
`DbsnpMyvariantAnnotator`, `ClinvarRsAnnotator`, `HgvsAnnotator`,
`SnpeffAnnotator`, `MafAnnotator`.

The many dbSNP annotators correspond to different ways of getting dbSNP data:
kind of 'scrapping' the web (not exactly, though), using Entrez service, and
using MyVariant.info service. Each gives a different piece of information, so
I included them all as separate entities.

All the annotators have the same API, just instantiate choosing a `cache` option
and call `annotate` with a list of identifiers.

## Cache

`anotamela` can use two types of cache: `RedisCache` and `PostgresCache`.

### Postgres
If you are to use Postgres, make sure you have it installed and create a 
database (e.g. `variants` would be an appropriate name) and a user with
privileges on it for `anotamela` to use. Then create a file (by default it will
be looked for in `~/.postgres_credentials.yml`) with your credentials:

```yaml
user: <your psql user>
pass: <your psql pass>
host: <where the server is, localhost if it's local>
port: <default for psql is 5432>
db: <the name of the db you created>
```

Then, you will have to `pip install sqlalchemy pyscopg2` and you will probably
have to `sudo apt-get install libpq-dev` for the last package. After this, you
are good to go.

```python
creds = '~/.psql_creds.yml'
annotator = DbsnpMyvariantAnnotator(cache='postgres', credentials_file=creds)
```

### Redis
You will need to `pip install redis`.

If you have a redis server running with default settings, you needn't do
anything extra. `anotamela` will try to connect at `localhost`'s port 6379.
If the server is running with non-default settings, specify those when
initializating each annotator:

```python
annotator = DbsnpWebAnnotator(cache='redis', host='localhost', port=5678)
```

### Sharing Cache between annotators

Several annotators can share the same Cache instance:

```python
psql_cache = PostgresCache(credentials_file='~/.psql_creds.yml')
dbsnp_ann = DbsnpMyvariantAnnotator(cache=psql_cache)
omim_ann = OmimVariantAnnotator(cache=psql_cache)
```

And you can later modify the cache of existing annotators:

```python
redis_cache = RedisCache(host='192.168.1.10')
dbsnp_ann.cache = redis_cache
omim_ann.cache = redis_cache
```

## Tests

`pytest -qx tests`

Tests take some time to run, since I didn't mock the web annotation. I
explicitely want to test the actual connection to each of the services and make
sure that the assumptions done in the code that parses the responses are still
sound.

