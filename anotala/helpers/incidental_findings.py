import re
from functools import lru_cache
from os.path import join, isfile
from datetime import datetime
from tempfile import gettempdir

import pandas as pd


@lru_cache()
def get_omim_incidental_genes_and_phenos():
    """
    Visits https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/ and gets the
    genes/phenos table as a pandas.DataFrame. Updates it once a month.
    """
    date_str = datetime.now().strftime('%Y_%m')
    fn = join(gettempdir(), 'ACMG_incidental_genes_{}.csv'.format(date_str))

    if not isfile(fn):
        tables = pd.read_html('https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/')
        assert len(tables) == 1  # Make sure there's only one table in the page
        df = tables[0]

        df.columns = df.iloc[0]
        df = df.drop(0).reset_index(drop=True)

        name = re.compile(r'(.+) ?\(')
        mim_id = re.compile(r'.+ ?\(MIM (\d+)\)')

        mapping = {'Disease name and MIM number': 'phenotype',
                   'Gene via GTR': 'gene'}
        for old_field, new_field in mapping.items():
            df[new_field] = df[old_field].str.extract(name, expand=False)
            df[new_field + '_MIM'] = df[old_field].str.extract(mim_id,
                                                               expand=False)
            df = df.drop(old_field, axis=1)

        no_gene = df['gene'].isnull()
        df.loc[no_gene, 'gene'] = \
            df.loc[no_gene, 'MedGen'].str.extract(r'MedGen\s+(\w+) \(', expand=False)
        df.loc[no_gene, 'gene_MIM'] = \
            df.loc[no_gene, 'MedGen'].str.extract(r'\(MIM (\d+)\)', expand=False)

        df = df.drop(['MedGen', 'Variations that maybe pathogenic'], axis=1)

        # Fix some missing values because colspan > 1
        last_row = None
        for i, row in df.iterrows():
            if any(row.isnull()):
                fixed = {}
                for key, value in row.fillna('').items():
                    if re.search(r'[A-Z0-9]+ \(MIM \d+\)$', value):
                        fixed['Gene via GTR'] = value
                    elif re.search(r'[a-z]+', value) and value not in ['ClinVar', 'MedGen']:
                        fixed['Disease name and MIM number'] = value
                    else:
                        fixed[key] = last_row[key]

                for fixed_key, fixed_value in fixed.items():
                    df.loc[i, fixed_key] = fixed_value

            last_row = row

        df.to_csv(fn, index=False)

    return pd.read_csv(fn, dtype={'gene_MIM': str, 'phenotype_MIM': str})


def is_incidental_gene(gene_mim_id):
    """
    Given a MIM ID for a gene entry, returns True if it should be reported
    as an incidental finding. Taken from:
    https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/
    """
    df = get_omim_incidental_genes_and_phenos()
    return str(gene_mim_id) in df['gene_MIM'].values


def is_incidental_pheno(pheno_mim_id):
    """
    Given a MIM ID for a gene entry, returns True if it should be reported
    as an incidental finding. Taken from:
    https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/
    """
    df = get_omim_incidental_genes_and_phenos()
    return str(pheno_mim_id) in df['phenotype_MIM'].values

