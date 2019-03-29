## Usages

### GRCh38 and GRCh37 positions for a list of rs IDs

```python
from anotala import DbsnpWebAnnotator

# You can increase these values if you're using Tor:
DbsnpWebAnnotator.BATCH_SIZE = 30
DbsnpWebAnnotator.SLEEP_TIME = 2

proxies = {'http': 'socks5://some-server:9050'}

dbsnp_web = DbsnpWebAnnotator('mysql', proxies=proxies)
annotations = dbsnp_web.annotate(rsids)

annotations_web
# => {
#   'GRCh37.p13_chrom': '7',
#   'GRCh37.p13_reverse': True,
#   'GRCh37.p13_start': '87230193',
#   'GRCh37.p13_stop': '87230193',
#   'GRCh38.p7_chrom': '7',
#   'GRCh38.p7_reverse': True,
#   'GRCh38.p7_start': '87600877',
#   'GRCh38.p7_stop': '87600877',
#   ...
# }
```


### PubMed PMIDs quick annotation from the command line

Given a list of PMIDs in a file (for instance, `/tmp/pubmeds.list`),
get profuse annotations in JSON format for each of them:

```bash
python anotala/cli.py pubmed annotate-pmids --pmids-file /tmp/pubmeds.list --cache mysql -o json
```

### Clinvar variants for a list of phenotypes, filtering by clinical significances

```python
from anotala.recipes import search_clinvar_variants_for_phenotypes

phenos = ['Diabetes melllitus', 'Alzheimer']

variants = search_clinvar_variants_for_phenotypes(
    pheno_terms=phenos, clinsigs=['pathogenic'], cache='mysql',
)

# You will get a list of ClinVar variants, each a dictionary full of info
# including dbSNP IDs, genes, associated phenotypes, HGVS names, etc.
# Basically, everything you see on the ClinVar web when you visit the
# details page for a variation.
```

You can do this from the command line also (WIP):

```bash
python anotala/cli.py clinvar search --phenos "Diabetes mellitus,Alzheimer"
--clinsigs pathogenic --cache mysql -o json
```

### Clinvar variants for a list of genes

```python
from anotala.recipes import clinvar_variants_for_genes

genes = ['BRCA2', 'CDH1', 'CDK4']
variant_ids, annotations = clinvar_variants_for_genes(genes, cache='mysql')

# You will get a list of ClinVar variant ids for those genes
# And a dictionary of annotations for each of those variants

# Optionally, you can get a DataFrame instead
df = clinvar_variants_for_genes(genes, cache='mysql', as_dataframe=True)
```

### GWAS Catalog list of traits for a given list of RS IDs

```python
from anotala.recipes import annotate_gwas_traits

proxies = {'http': 'socks5://localhost:9050'} # TOR instance!

rsids = ['rs123', 'rs234']

annotate_gwas_traits(rsids, cache='dict', proxies=proxies)
# => {'rs123': ['Eye color traits'], 'rs268': ['Metabolic syndrome']}
```

### Annotation pipeline

My most common pattern involves using the whole annotation pipeline as boxed
by the package. I usually do this in Jupyter, and making use of other
package called [`project`](http://github.com/biocodices/project):

```python
WORKDIR = 'some_descriptive_name_for_the_project'

from project import Project
from anotala import AnnotationPipeline

pj = Project(WORKDIR)

# Typical place of a local Tor running:
proxies = {'http': 'socks5://localhost:9050'}

# The result of your variant calling:
genotypes = pj.results_file('genotypes.vcf.gz')

# I use MySQL as a cache for the web annotations.
# This will look for host/user/pass/db in a ~/.mysql_credentials.yml file.
annotation_pipe = AnnotationPipeline(cache='mysql', proxies=proxies)

# The annotations itself will take some time:
annotation_pipe.run_from_vcf(genotypes)

# => variants are stored in annotation_pipe.rs_variants
# => genes are stored in annotation_pipe.gene_annotations

# Dump the results
ann_variants_1kg_json = \
    pj.dump_df_as_json(annotation_pipe.rs_variants, 'annotated_rs_variants')
ann_genes_1kg_json = \
    pj.dump_df_as_json(annotation_pipe.gene_annotations, 'annotated_genes')
```

This produces two pandas DataFrames that I dump as JSON files to be used
downstream in a reports generation pipeline.

However, the annotators can be used separately, like this:

```python
from anotala import DbsnpMyvariantAnnotator, OmimVariantAnnotator

dbsnp_annotator = DbsnpMyvariantAnnotator(cache='redis')
dbsnp_annotator.annotate(['rs268', 'rs123'])
# => {'rs268': { ... }, 'rs123': { ... }}

omim_annotator = OmimVariantAnnotator(cache='postgres')
omim_annotator.annotate(['605557.0001', '605557.0003'])
# => {'605557.0001': { ... }, '605557.003': { ... }}

# Just fetch from cache, don't use the Internets:
omim_annotator.annotate(rs_list, use_web=False)

# Fetch everything from the Net, don't use the cache.
# This will also update the existing cache entries:
omim_annotator.annotate(rs_list, use_cache=False)
```

There are plenty of annotators: `DbsnpWebAnnotator`, `DbsnpEntrezAnnotator`,
`DbsnpMyvariantAnnotator`, `ClinvarRsAnnotator`, `HgvsAnnotator`,
`SnpeffAnnotator`, `MafAnnotator`, `OmimGeneAnnotator`, `OmimVariantAnnotator`,
`GeneEntrezAnnotator`, `MygeneAnnotator`, `PubmedAnnotator`, `UniprotAnnotator`.

As a general rule, SNP annotators take rs IDs, gene annotators take gene Entrez
IDs, and OMIM annotators take MIM IDs.

The many dbSNP annotators correspond to different ways of getting dbSNP data:
kind of 'scrapping' the web, using Entrez service, and using MyVariant.info
service. Each gives a different piece of information, so I included them as
separate entities.

All the annotators have the same API, just instantiate choosing a `cache` option
and call `annotate` with a list of identifiers.

A common pipeline would involve annotating a list of rs IDs with DbSNP and
ClinVar, then taking the OMIM IDs provided by ClinVar and using those to
annotate with OMIM.

## Cache

`anotala` can use different types of cache: `RedisCache`, `MysqlCache`,
`PostgresCache`, `DictCache`.

### MySQL

Make sure you have mysql installed, create a database (e.g. `anotala_cache`)
and a user with privileges on that database. Then create a YAML file (by default
it will be looked for in `~/.mysql_credentials.yml`) with the credentials:

```yaml
host: <where the server is, "localhost" if it's local>
user: <your mysql user, e.g. "juan">
pass: <your mysql pass, e.g. "mypass">
port: <usually 3306 for mysql>
db: <the name of an *existing* db, e.g. "anotala_cache">
driver: <for instance, "mysql+pymysql", if you have this one installed>
```

### Postgres

If you are to use Postgres, make sure you have it installed and create a 
database (e.g. `anotala_cache` would be an appropriate name) and a user with
privileges on it for `anotala` to use. Then create a file (by default it will
be looked for in `~/.postgres_credentials.yml`) with your credentials:

```yaml
user: <your psql user>
pass: <your psql pass>
host: <where the server is, 'localhost' if it's local>
port: <default for psql is 5432>
db: <the name of an *existing* db, e.g. 'variants'>
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
anything extra. `anotala` will try to connect at `localhost`'s port 6379.
If the server is running with non-default settings, specify those when
initializating each annotator:

```python
annotator = DbsnpWebAnnotator(cache='redis', host='localhost', port=5678)
```

### Dict

The `DictCache` (named `dict`) is meant for testing or for a quick use
when you don't have any database software installed, like MySQL, Redis, etc.
With this option, annotators will use a Python dictionary as their "cache",
which means there is no persistance of any kind: all fetched data is lost
whenever the Python session is ended and is not even shared between annotators.

Not the best option, but it comes handy once in a while, in some settings.

### Sharing Cache between annotators

Several annotators can share the same Cache instance:

```python
psql_cache = PostgresCache(credentials_file='~/.psql_creds.yml')
dbsnp_ann = DbsnpMyvariantAnnotator(cache=psql_cache)
omim_ann = OmimVariantAnnotator(cache=psql_cache)
```

And you can later modify the cache of existing annotators:

```python
dbsnp_ann = DbsnpMyvariantAnnotator(cache='postgres')
dbsnp_ann.annotate('rs268')  # Data will be written to PostgreSQL

redis_cache = RedisCache()
dbsnp_ann.cache = redis_cache
dbsnp_ann.annotate('rs268')  # Data will be written to Redis now
```

## Tests

`pytest -qx tests`

Tests take some time to run, since I didn't mock the web annotation. I
explicitely want to test the actual connection to each of the services and make
sure that the assumptions done in the code that parses the responses are still
sound.

Tests are split in two directories according to wether they use the web
(i.e. they are slow to run) or not, so you can run the non-web tests quickly:
`pytest tests`, and later just the web tests: `pytest tests_that_use_web`.

