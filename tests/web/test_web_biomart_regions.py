from anotamela.annotators import BiomartRegionsAnnotator


def test_biomart_regions_from_web(proxies):
    annotator = BiomartRegionsAnnotator(cache='dict', proxies=proxies)

    region = '12:21953990:21953990:1'
    snp = {'chrom_end_g37': 21953990,
           'chrom_g37': '12',
           'chrom_start_g37': 21953990,
           'chrom_strand_g37': 1,
           'rsid': 'rs143310355',
           'source': 'dbSNP'}

    result = annotator.annotate_one(region)
    assert result == [snp]

