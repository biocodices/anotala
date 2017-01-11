import re

from anotamela.annotators.base_classes import MyVariantAnnotator


SNP_RE = re.compile(r'(?P<old_allele>[ATCG])>?(?P<new_allele>[A-Z|=])')


class MafAnnotator(MyVariantAnnotator):
    SOURCE_NAME = 'maf'
    SCOPES = 'dbsnp.rsid dbnsfp.rsid'.split()
    FIELDS = ('dbsnp.gmaf dbnsfp.1000gp3 dbnsfp.exac dbnsfp.esp6500 '
              'dbnsfp.twinsuk.af cadd.1000g').split()

    @staticmethod
    def _parse_hit(hit):
        hgvs = hit['_id']
        if SNP_RE.search(hgvs):
            allele = SNP_RE.search(hgvs).group('new_allele')
        else:
            allele = hgvs

        dicts = (
            ('1000gp3', hit.get('dbsnp', {})),
            ('1000gp3', hit.get('dbnsfp', {}).get('1000gp3', {})),
            ('esp6500', hit.get('dbnsfp', {}).get('esp6500', {})),
            ('twinsuk', hit.get('dbnsfp', {}).get('twinsuk', {})),
            ('exac', hit.get('dbnsfp', {}).get('exac', {})),
            ('cadd', hit.get('cadd', {}).get('1000g', {}))
        )

        annotation = {}
        for source_name, dic in dicts:
            for key, value in dic.items():
                compound_key = '{}_{}_{}'.format(source_name, key, allele)
                # Doesn't make sense to have a MAF with 10 decimals, round it:
                annotation[compound_key] = round(value, 4)

        # Remove allele counts, leave frequency
        return {re.sub('_af\b', '', k): v for k, v in annotation.items()
                if not k.endswith('_ac')}

