from anotamela import DbsnpWebAnnotator
import pytest


@pytest.fixture
def raw_annotation():
    return {
        'snp_id': '1',
        'assembly': {
            'GRCh37.p13': [
                {
                    'groupTerm': 'Primary_Assembly', 'contigLabel': '',
                    'snp2chrOrien': '1',
                    'chr': '1',
                    'chrPosFrom': '100',
                    'chrPosTo': '200',
                    'geneModel': [
                        {
                            'geneId': '2',
                            'geneSymbol': 'GENE-2',
                            'variation_allele': [
                                {'fxnName': 'UTR-5'}
                            ]
                        },
                    ]
                 },  # Should keep this one

                # Ignore PAR (Pseudo-autosomal region)
                {'groupTerm': 'Primary_Assembly',
                 'contigLabel': 'PAR',
                 'snp2chrOrien': '1'},

                # Ignore non Primary-Assembly
                {'groupTerm': 'PATCHES',
                 'contigLabel': '',
                 'snp2chrOrien': '0'},
            ],

            'GRCh38.p7': [
                {
                    'groupTerm': 'Primary_Assembly',
                    'contigLabel': '',
                    'snp2chrOrien': '0',
                    'chr': '1',
                    'chrPosFrom': '110',
                    'chrPosTo': '210',
                    'geneModel': [
                        {'geneId': '1', 'geneSymbol': 'GENE-1'},
                        {'geneId': '1', 'geneSymbol': 'GENE-1'},

                        {
                            'codonPos': '1',
                            'contigAcc': 'NT_167190.2',
                            'contig_allele': {
                                'aaAbbrev': 'Tyr',
                                'aaCode': 'Y',
                                'allele': 'T',
                                #  'codon': 'TAC',
                                #  'fxnClass': '8',
                                #  'fxnName': 'cds-reference'
                            },
                            'geneId': '2',
                            'geneSymbol': 'GENE-2',
                            'proteinAcc': 'NP_001157289.1',
                            'proteinPos': '408',
                            'snp2mRNAOrien': '1',
                            'variation_allele': [
                                {
                                    'aaAbbrev': 'His',
                                    'aaCode': 'H',
                                    'allele': 'C',
                                    'codon': 'CAC',
                                    'fxnName': 'missense'
                                }
                            ],
                        },
                    ],
                },

                {
                    'groupTerm': 'PATCHES',
                    'contigLabel': '',
                    'snp2chrOrien': '1'
                },
            ],
        }
    }


def test_parse_annotation(raw_annotation):

    parsed = DbsnpWebAnnotator._parse_annotation(raw_annotation)

    assert parsed['consequences_per_gene'] == {'GENE-2': ['UTR-5', 'missense']}
    # Explicitely test booleans, not truthiness:
    assert parsed['GRCh37.p13_reverse'] is True
    assert parsed['GRCh37.p13_chrom'] == '1'
    assert parsed['GRCh37.p13_start'] == 100
    assert parsed['GRCh37.p13_stop'] == 200

    assert parsed['GRCh38.p7_chrom'] == '1'
    assert parsed['GRCh38.p7_start'] == 110
    assert parsed['GRCh38.p7_stop'] == 210
    assert parsed['GRCh38.p7_reverse'] is False
    assert sorted(parsed['GRCh38.p7_gene_symbols']) == ['GENE-1', 'GENE-2']
    assert sorted(parsed['GRCh38.p7_gene_entrez_ids']) == [1, 2]

    assert parsed['rs_id'] == 'rs1'

    incomplete_annotation = {
        'snp_id': '2',
        'assembly': {
            'GRCh37.p13': [{'groupTerm': 'ALT_LOCI'}],
        }
    }

    # Shouldn't break:
    parsed = DbsnpWebAnnotator._parse_annotation(incomplete_annotation)

    assert parsed['rs_id'] == 'rs2'


def test_extract_consequences_per_gene_from_gene_models():
    gene_models = [
        {
            'geneSymbol': 'GENE-1',
            'variation_allele': [
                {'fxnName': 'MISSENSE'},
                {'fxnName': 'UTR-5'}
            ]
        },
        {
            'geneSymbol': 'GENE-2',
            'variation_allele': [
                {'fxnName': 'SYNONYMOUS'},
            ]
        },
    ]
    result = DbsnpWebAnnotator._extract_consequences_per_gene_from_gene_models(gene_models)
    assert result == {
        'GENE-1': ['MISSENSE', 'UTR-5'],
        'GENE-2': ['SYNONYMOUS'],
    }
