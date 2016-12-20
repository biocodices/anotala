from os.path import expanduser, isfile
from functools import lru_cache

import requests
import pandas as pd


@lru_cache()
def mim_to_gene_df():
    """
    Returns a DataFrame with:
        - MIM ID
        - MIM Entry Type
        - Entrez Gene ID
        - HGNC Symbol (Gene Symbol)
        - Ensembl Gene ID
    """
    cache_file = expanduser('~/.mim2gene.txt')

    if not isfile(cache_file):
        url = 'https://omim.org/static/omim/data/mim2gene.txt'
        response = requests.get(url)
        if response.ok:
            with open(cache_file, 'w') as f:
                f.write(response.text)
        else:
            response.raise_for_status()

    fields = 'mim_id entry_type entrez_id gene_symbol ensembl_ids'.split()
    df = pd.read_table(cache_file, comment='#', names=fields,
                       dtype={'entrez_id': str, 'mim_id': str})

    return df.set_index('mim_id')


@lru_cache()
def mim_to_gene(id_=None):
    """Maps OMIM IDs to Entrez gene IDs. Pass and ID to get the MIM or no args
    to get the whole dictionary."""
    df = mim_to_gene_df()
    dic = {k: v for k, v in df['entrez_id'].dropna().items()}
    return dic[id_] if id_ else dic


@lru_cache()
def gene_to_mim(id_=None):
    """Maps Entrez gene IDs to OMIM IDs. Pass and ID to get the MIM or no args
    to get the whole dictionary."""
    df = mim_to_gene_df()
    dic = {v: k for k, v in df['entrez_id'].dropna().items()}
    return dic[id_] if id_ else dic

