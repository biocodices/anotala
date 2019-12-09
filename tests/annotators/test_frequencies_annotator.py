from anotala import FrequenciesAnnotator


freq, count = 0.1234010101, 100  # Fake frequency and allele count data
rounded_freq = 0.123401  # Rounded!


def test_frequencies_annotator_parse_hit():
    hit = {
        'allele': 'G',  # This datum comes from a previous parsing

        'cadd': {'1000g': {'af': freq,
                           'afr': freq,
                           'amr': freq,
                           'eur': freq}},
        'gnomad_exome': {'af': {'af_afr': freq}},
        'gnomad_genome': {'af': {'af_afr': freq}},

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
        'African (AFR)',
        'American (AMR)',
        'European (EUR)',
    ]
    for population in cadd_populations:
        assert G['CADD_1000g'][population] == rounded_freq

    esp_populations = [
        'African American (AA)',
        'European American (EA)',
    ]
    for population in esp_populations:
        assert G['dbNSFP_ESP6500'][population] == rounded_freq

    kg3_populations = [
        'General',
        'African (AFR)',
        'American (AMR)',
        'East Asian (EAS)',
        'European (EUR)',
        'South Asian (SAS)',
    ]
    for population in kg3_populations:
        assert G['dbNSFP_1000gp3'][population] == rounded_freq

    exac_populations = [
        'General (ADJ)',
        'General',
        'African (AFR)',
        'American (AMR)',
        'East Asian (EAS)',
        'Finnish (FIN)',
        'Non-Finnish European (NFE)',
        'South Asian (SAS)',
    ]
    for population in exac_populations:
        assert G['dbNSFP_ExAC'][population] == rounded_freq

    assert G['dbNSFP_twinsUK']['General'] == rounded_freq

    assert G['gnomad_exome']['African (AFR)'] == rounded_freq
    assert G['gnomad_genome']['African (AFR)'] == rounded_freq


def test_parse_hit_with_few_data():
    hit = {
        'allele': 'G',
        'cadd': {'1000g': {'afr_af': freq}}
    }

    r = FrequenciesAnnotator._parse_hit(hit)

    assert r == {'G': {'CADD_1000g': {'African (AFR)': rounded_freq}}}


def test_parse_annotations_hook():
    annotations = [
        {'A': {'foo': 'bar'}, 'G': {'baz': 'qux'}},
        {'G': {'spam': 'eggs'}},
    ]

    result = FrequenciesAnnotator._parse_annotations_hook(annotations)
    assert result['G'] == {'baz': 'qux', 'spam': 'eggs'}
    assert result['A'] == {'foo': 'bar'}

