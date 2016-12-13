from anotamela.annotators.base_classes import MyVariantAnnotator


class MafAnnotator(MyVariantAnnotator):
    SOURCE_NAME = 'maf'
    SCOPES = 'dbsnp.rsid dbnsfp.rsid'.split()
    FIELDS = 'dbsnp.gmaf dbnsfp.1000gp3 dbnsfp.exac dbnsfp.esp6500'.split()

    @staticmethod
    def _parse_hit(hit):
        dicts = (
            ('1000gp3', hit.get('dbsnp', {})),
            ('1000gp3', hit.get('dbnsfp', {}).get('1000gp3', {})),
            ('esp6500', hit.get('dbnsfp', {}).get('esp6500', {})),
            ('exac', hit.get('dbnsfp', {}).get('exac', {}))
        )

        annotation = {}
        for source_name, dic in dicts:
            for key, value in dic.items():
                compound_key = '{}_{}'.format(source_name, key)
                # Doesn't make sense to have a MAF with 10 decimals, round it:
                annotation[compound_key] = round(value, 4)

        # Remove allele counts, leave frequency
        return {k: v for k, v in annotation.items()
                if not k.endswith('_ac')}
