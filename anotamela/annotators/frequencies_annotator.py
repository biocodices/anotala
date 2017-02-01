from collections import defaultdict

from anotamela.annotators.base_classes import MyVariantAnnotator


class FrequenciesAnnotator(MyVariantAnnotator):
    SOURCE_NAME = 'frequencies'
    SCOPES = 'dbsnp.rsid dbnsfp.rsid'.split()
    FIELDS = ('dbsnp.alleles dbnsfp.1000gp3 dbnsfp.exac '
              'dbnsfp.esp6500 dbnsfp.twinsuk.af cadd.1000g').split()

    @classmethod
    def _parse_hit(cls, hit):
        allele = hit['allele']
        frequencies = defaultdict(dict)

        for entry in hit.get('dbsnp', {}).get('alleles', []):
            if 'freq' in entry:
                freq = cls._parse_frequency(entry['freq'])
                frequencies[entry['allele']] = \
                    {'dbSNP': {'General': freq}}

        sources = {
            'dbNSFP_1000gp3': hit.get('dbnsfp', {}).get('1000gp3', {}),
            'dbNSFP_ESP6500': hit.get('dbnsfp', {}).get('esp6500', {}),
            'dbNSFP_twinsUK': hit.get('dbnsfp', {}).get('twinsuk', {}),
            'dbNSFP_ExAC': hit.get('dbnsfp', {}).get('exac', {}),
            'CADD_1000g': hit.get('cadd', {}).get('1000g', {}),
        }

        for source_name, source_dict in sources.items():
            frequencies[allele][source_name] = info = {}
            for population, freq in source_dict.items():
                population = cls._parse_population(population)
                if population:
                    info[population] = cls._parse_frequency(freq)

        return dict(frequencies)

    @staticmethod
    def _parse_frequency(freq):
        """Round the allele frequency."""
        return round(freq, 4)

    @staticmethod
    def _parse_population(population):
        """
        Translate tags like 'afr_af' to a description like 'African'.
        Return None for allele count tags like 'afr_ac'.
        """
        if population == 'ac' or population.endswith('_ac'):
            return  # Ignore allele counts (ac) data

        population = population.replace('_af', '')
        abbreviations = {
            'af': 'General',
            'ea': 'European American',
            'aa': 'African American',
            'adj': 'General (Adj.)',
            'amr': 'American (native)',
            'eur': 'European',
            'afr': 'African',
            'eas': 'East Asian',
            'sas': 'South Asian',
            'fin': 'Finnish',
            'nfe': 'Non-Finnish European',
            'oth': 'Other',
        }
        return abbreviations.get(population, population.upper())

