## Annotation

## Usage
```python
from anotamela import DbsnpWebAnnotator

dbsnp_web = DbsnpWebAnnotator(cache='redis')  # or cache='postgres'
dbsnp_web.annotate(['rs268', 'rs123'])
# => {'rs268': { ... }, 'rs123': { ... }}

dbsnp_entrez = DbsnpEntrezAnnotator(cache='postgres')
dbsnp_entrez.annotate(['rs268', 'rs123'])
# => {'rs268': { ... }, 'rs123': { ... }}

# You can tune the parallelization options
# Make 5 parallel requests each time, sleep 10s between batches:
dbsnp_web.annotate(rs_list, parallel=5, sleep_time=10)

# Just fetch from cache, no web:
dbsnp_web.annotate(rs_list, use_web=False)

# Fetch everything from web, don't use the cache.
# This will update the existing cache entries also:
dbsnp_web.annotate(rs_list, use_cache=False)
```

There are plenty of annotators: `DbsnpWebAnnotator`, `DbsnpEntrezAnnotator`,
`ClinvarRsAnnotator`, `HgvsAnnotator`, `SnpeffAnnotator`.
All of them have the same API, just instantiate choosing a `cache` option and
call `annotate` with a list of identifiers.

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

### Redis
You will need to `pip install redis`.

If you have a redis server running with default settings, you needn't do
anything extra. `anotamela` will try to connect at `localhost`'s port 6379.
If the server is running with non-default settings, specify those when
initializating each annotator:

```python
DbsnpWebAnnotator(cache='redis', host='localhost', port=5678)
```

## Tests

`pytest -qx tests`

Tests take some time to run, since I didn't mock the web annotation. I
explicitely want to test the actual connection to each of the services and make
sure that the assumptions done in the code are still sound when testing.

