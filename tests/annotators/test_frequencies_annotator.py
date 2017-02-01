from anotamela import FrequenciesAnnotator


def test_maf_annotator_parse_hit():
    freq, count = 0.123401, 100  # Fake frequency and allele count data
    rounded_freq = 0.1234  # Rounded!

    hit = {
        'allele': 'G',  # This datum comes from a previous parsing

        'cadd': {'1000g': {'af': freq,
                           'afr': freq,
                           'amr': freq,
                           'eur': freq}},

        'dbnsfp': {'1000gp3': {'ac': count,
                               'af': freq,
                               'afr_ac': count,
                               'afr_af': freq,
                               'amr_ac': count,
                               'amr_af': freq,
                               'eas_ac': count,
                               'eas_af': freq,
                               'eur_ac': count,
                               'eur_af': freq,
                               'sas_ac': count,
                               'sas_af': freq},

                   'esp6500': {'aa_ac': count,
                               'aa_af': freq,
                               'ea_ac': count,
                               'ea_af': freq},

                   'exac': {'ac': count,
                            'adj_ac': count,
                            'adj_af': freq,
                            'af': freq,
                            'afr_ac': count,
                            'afr_af': freq,
                            'amr_ac': count,
                            'amr_af': freq,
                            'eas_ac': count,
                            'eas_af': freq,
                            'fin_ac': count,
                            'fin_af': freq,
                            'nfe_ac': count,
                            'nfe_af': freq,
                            'sas_ac': count,
                            'sas_af': freq},

                   'twinsuk': {'af': freq}},

        'dbsnp': {'alleles': [{'allele': 'A', 'freq': freq},
                              {'allele': 'G', 'freq': freq},
                              {'allele': 'T'}],  # has no 'freq'
                  'gmaf': freq}
    }

    # Produces a per-allele dictionary with dictionaries of
    # population->frequency per source. Rounds the floats with 4 precision:

    r = FrequenciesAnnotator._parse_hit(hit)

    assert 'G' in r
    assert 'A' in r
    assert 'T' not in r

    assert r['A']['dbSNP']['General'] == rounded_freq

    G = r['G']

    assert G['dbSNP']['General'] == rounded_freq

    cadd_populations = [
        'General',
        'African',
        'American (native)',
        'European',
    ]
    for population in cadd_populations:
        assert G['CADD_1000g'][population] == rounded_freq

    esp_populations = [
        'African American',
        'European American',
    ]
    for population in esp_populations:
        assert G['dbNSFP_ESP6500'][population] == rounded_freq

    kg3_populations = [
        'General',
        'African',
        'American (native)',
        'East Asian',
        'European',
        'South Asian',
    ]
    for population in kg3_populations:
        assert G['dbNSFP_1000gp3'][population] == rounded_freq

    exac_populations = [
        'General (Adj.)',
        'General',
        'African',
        'American (native)',
        'East Asian',
        'Finnish',
        'Non-Finnish European',
        'South Asian',
    ]
    for population in exac_populations:
        assert G['dbNSFP_ExAC'][population] == rounded_freq

    assert G['dbNSFP_twinsUK']['General'] == rounded_freq

