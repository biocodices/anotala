from anotamela.annotators import AnnotatorWithCache, MyVariantAnnotator


class MafAnnotator(MyVariantAnnotator, AnnotatorWithCache):
    SOURCE_NAME = 'maf'
    SCOPES = 'dbnsfp.rsid'
    FIELDS = 'dbsnp.gmaf dbnsfp.1000gp3 dbnsfp.exac dbnsfp.esp6500'.split()

    @staticmethod
    def _parse_annotation(raw_annotation):
        dicts = (
            ('1000gp3', raw_annotation.get('dbsnp', {})),
            ('1000gp3', raw_annotation.get('dbnsfp', {}).get('1000gp3', {})),
            ('esp6500', raw_annotation.get('dbnsfp', {}).get('esp6500', {})),
            ('exac', raw_annotation.get('dbnsfp', {}).get('exac', {}))
        )

        annotation = {}
        for source_name, dic in dicts:
            for key, value in dic.items():
                compound_key = '{}_{}'.format(source_name, key)
                # Doesn't make sense to have a MAF with 10 decimals, round it:
                annotation[compound_key] = round(value, 4)

        return {k: v for k, v in annotation.items()
                if not k.endswith('_ac')}  # Remove allele counts, leaves freq

