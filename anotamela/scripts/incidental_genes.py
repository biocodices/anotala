#!/usr/bin/env python

import re
from os.path import dirname, basename, join

import pandas as pd


url = 'https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/'
df = pd.read_html(url)[0]
df.columns = df.iloc[0]
df = df.drop(0).reset_index(drop=True)

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

# Organize the dataframe
mapping = {'Disease name and MIM number': 'phenotype',
           'Gene via GTR': 'gene'}

for old_field, new_field in mapping.items():
    df[new_field] = df[old_field].str.extract(r'(.+) \(', expand=False)
    df[new_field + '_MIM'] = df[old_field].str.extract(r'.+ \(MIM (\d+)\)',
                                                       expand=False)
    df.drop(old_field, axis=1, inplace=True)

df = df.drop(['MedGen', 'Variations that may be pathogenic'], axis=1)

target = join(dirname(__file__), basename(__file__).replace('.py', '.csv'))
df.to_csv(target, index=None)
print('CSV with incidental genes written to: "{}"'.format(target))

