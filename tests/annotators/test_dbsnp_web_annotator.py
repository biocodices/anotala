from anotamela import DbsnpWebAnnotator


def test_parse_annotation():
    annotation = {
        'assembly': {
            'GRCh37.p13': [
                {
                    'groupTerm': 'Primary_Assembly', 'contigLabel': '',
                    'snp2chrOrien': '1',
                    'chr': '1',
                    'chrPosFrom': '100',
                    'chrPosTo': '200',
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
                    'groupTerm': 'Primary_Assembly', 'contigLabel': '',
                    'snp2chrOrien': '0',
                    'chr': '1',
                    'chrPosFrom': '110',
                    'chrPosTo': '210',
                    'geneModel': [
                        {'geneSymbol': 'GENE-1'},
                        {'geneSymbol': 'GENE-1'},
                    ]
                 },

                {'groupTerm': 'PATCHES', 'contigLabel': '',
                 'snp2chrOrien': '1'},
            ],
        }
    }

    parsed = DbsnpWebAnnotator._parse_annotation(annotation)

    # Explicitely test booleans, not truthiness:
    assert parsed['GRCh37.p13_reverse'] is True
    assert parsed['GRCh37.p13_chrom'] == '1'
    assert parsed['GRCh37.p13_start'] == 100
    assert parsed['GRCh37.p13_stop'] == 200

    assert parsed['GRCh38.p7_chrom'] == '1'
    assert parsed['GRCh38.p7_start'] == 110
    assert parsed['GRCh38.p7_stop'] == 210
    assert parsed['GRCh38.p7_reverse'] is False
    assert parsed['GRCh38.p7_gene_symbols'] == ['GENE-1']

    annotation = {
        'assembly': {
            'GRCh37.p13': [{'groupTerm': 'ALT_LOCI'}],
        }
    }

    # Shouldn't break:
    DbsnpWebAnnotator._parse_annotation(annotation)
