from anotamela.helpers import ClinvarVCFParser

import pandas as pd
import pytest


def path_to_vcf():
    return pytest.helpers.file('clinvar.vcf')


@pytest.fixture
def parser():
    return ClinvarVCFParser()


def test_read_file(parser):
    result = parser.read_file(path_to_vcf())
    assert len(result) == 3


def test_parse_data(parser):
    input_df = pd.DataFrame([
        {'id': 'Var-1',
         'info': {'RS': '1',
                  'CLNVC': 'Type-1',
                  'CLNHGVS': 'name-1',
                  'CLNDISDB': 'OMIM:OMIM0|OMIM:OMIM1,Human_Phenotype_Ontology:HP:1,UnknownDB:1',
                  'CLNSIG': 'Sig-1',
                  'CLNVI': 'Section-1,Center-1:CS1,OMIM_Allelic_Variant:100.001,PharmGKB:P1,UniProtKB_(protein):P1#VAR_1',
                  'GENEINFO': 'GENE-1:1|GENE-2:2',
                  'CLNVCSO': 'SO:1',
                  'ORIGIN': '1',
                  'MC': 'SO:0|consequence-1,SO:1|consequence-2',
                  'FOO': 'foo'}},
        {'info': {'RS': '2',
                  'ORIGIN': '2',
                  'CLNHGVS': 'name-2',
                  'BAR': 'bar'}},
        {'info': {'RS': '3',
                  'ORIGIN': 'UNKNOWN-KEY',
                  'CLNHGVS': 'name-3'}},
    ])
    df = parser.parse_data(input_df)

    # renames columns
    assert df.loc[0, 'variation_id'] == 'Var-1'
    assert df.loc[0, 'rs_id'] == 'rs1'
    assert df.loc[0, 'clinvar_hgvs'] == 'name-1'
    assert df.loc[0, 'source_info'] == input_df.loc[0, 'info']
    assert df.loc[0, 'clinical_significance'] == 'Sig-1'
    assert df.loc[0, 'variant_type'] == 'Type-1'
    assert df.loc[0, 'sequence_ontology_id'] == 'SO:1'

    # removes unnecessary columns
    assert 'RS' not in df
    assert 'CLNVI' not in df
    assert 'GENEINFO' not in df
    assert 'ORIGIN' not in df
    assert 'MC' not in df

    # adds "rs" to RS numbers
    assert df['rs_id'].values.tolist() == ['rs1', 'rs2', 'rs3']

    # extract keys as new columns
    assert df['FOO'].values.tolist() == ['foo', None, None]
    assert df['BAR'].values.tolist() == [None, 'bar', None]

    # extracts external db disease IDs
    assert df.loc[0, 'disease_ids'] == [
        'OMIM:OMIM0',
        'OMIM:OMIM1',
        'Human_Phenotype_Ontology:HP:1',
        'UnknownDB:1',
    ]
    assert df.loc[1, 'disease_ids'] == []

    # parses and extracts external sources IDs
    assert df.loc[0, 'sources'] == {
        'Section-1,Center-1': 'CS1',
        'OMIM_Allelic_Variant': '100.001',
        'PharmGKB': 'P1',
        'UniProtKB_(protein)': 'P1#VAR_1',
    }
    assert df.loc[0, 'omim_variant_id'] == '100.001'
    assert df.loc[0, 'pharmgkb_id'] == 'P1'
    assert df.loc[0, 'uniprot_variant_id'] == 'P1#VAR_1'

    assert df.loc[1, 'sources'] == None
    assert df.loc[1, 'omim_variant_id'] == None

    # parses gene info
    assert df.loc[0, 'genes'] == {'GENE-1': 1, 'GENE-2': 2}
    assert df.loc[0, 'gene_symbols'] == ['GENE-1', 'GENE-2']
    assert df.loc[0, 'gene_entrez_ids'] == [1, 2]

    assert df.loc[1, 'gene_symbols'] == []
    assert df.loc[1, 'gene_entrez_ids'] == []

    # parse origin
    assert df['origin'].values.tolist() == ['germline', 'somatic', 'UNKNOWN-KEY']

    # extract molecular consequences
    assert df.loc[0, 'molecular_consequences'] == ['consequence-1', 'consequence-2']
    assert df.loc[1, 'molecular_consequences'] == []

