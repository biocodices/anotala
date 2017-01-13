from anotamela import MafAnnotator


def test_maf_annotator_parse_hit():
    raw = {'_id': 'chr2:g.73679956C>G',
           '_score': 16.956917,
           'dbnsfp': {'exac': {'ac': 36,
                               'adj_ac': 35,
                               'adj_af': 0.0002962,
                               'af': 0.0002981,
                               'nfe_ac': 35,
                               'nfe_af': 0.0005379}},
           'dbsnp': {'gmaf': 0.01478},
           'query': 'rs28730854'}

    # Produces a 1-level dictionary with just the frequencies (not the counts),
    # Includes the allele in the keys,
    # Rounds the floats:
    assert MafAnnotator._parse_hit(raw) == {'1000gp3_gmaf_G': 0.0148,
                                            'exac_adj_af_G': 0.0003,
                                            'exac_af_G': 0.0003,
                                            'exac_nfe_af_G': 0.0005}

