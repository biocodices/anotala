from os.path import dirname, join

import pytest

from anotamela import Pipeline


TEST_PARAMS = ['test.vcf', 'test.vcf.gz']

@pytest.mark.parametrize('vcf_filename', TEST_PARAMS)
def test_pipeline(vcf_filename):
    # Test the pipeline using the web to annotate
    pipeline = Pipeline()
    proxies = {'http': 'socks5://beleriand.local:9150'}
    pipeline.run(vcf_path=_test_file(vcf_filename),
                 cache='mock_cache', use_cache=False, proxies=proxies)

    _test_pipeline_result(pipeline)

    # Test the pipeline again, now using the cache built in the test above
    new_pipeline = Pipeline()
    new_pipeline.run(vcf_path=_test_file(vcf_filename),
                     cache=pipeline.cache, use_web=False)

    _test_pipeline_result(pipeline)

def _test_pipeline_result(pipeline):
    assert 'rs28935490' in pipeline.rs_variants['id']

    vcf_fields = 'chrom pos id ref alt qual filter info'.split()
    assert 'esv3585040' in pipeline.other_variants['id']

    for vcf_field in vcf_fields:
        assert vcf_field in pipeline.other_variants.iloc[0]

    annotation_fields = ('clinvar_rs snpeff_myvariant maf hgvs '
                         'dbsnp_myvariant dbsnp_entrez entrez_gene_ids '
                         'omim_entries uniprot_entries').split()

    for row in pipeline.rs_variants.iterrows():
        for field in vcf_fields + annotation_fields:
            assert row[field]

    for gene_id in ['2717', '100529097']:
        assert gene_id in pipeline.gene_annotations['entrez_id']

def _test_file(filename):
    return join(dirname(__file__), 'files', filename)
