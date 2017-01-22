import pytest

from anotamela import AnnotationPipeline
from anotamela.cache import DictCache


TEST_PARAMS = [
    # file, test web?
    ('test.vcf', True),
    ('test.vcf.gz', False),
]

@pytest.fixture(scope='module')
def cache():
    return DictCache()


@pytest.mark.parametrize('vcf_filename,test_web', TEST_PARAMS)
def test_pipeline(proxies, vcf_filename, test_web, cache):
    # Test the pipeline using the web to annotate and build the cache
    if test_web:
        web_pipeline = AnnotationPipeline(cache=cache, use_cache=False,
                                          proxies=proxies)
        web_pipeline.run_from_vcf(vcf_path=pytest.helpers.file(vcf_filename))
        _test_pipeline_result(web_pipeline)

    # Test the pipeline again, now using the cache built previously
    cache_pipeline = AnnotationPipeline(cache=cache, use_web=False)
    cache_pipeline.run_from_vcf(vcf_path=pytest.helpers.file(vcf_filename))
    _test_pipeline_result(cache_pipeline)


def _test_pipeline_result(pipeline):
    assert 'esv3585040' in pipeline.other_variants['id'].values

    rsid = 'rs28935490'
    assert rsid in pipeline.rs_variants['id'].values

    vcf_fields = 'chrom pos id ref alt qual filter info'.split()
    annotation_fields = ('clinvar_entries snpeff_myvariant maf hgvs dbnsfp '
                         'dbsnp_myvariant entrez_gene_ids dbsnp_entrez '
                         'omim_entries uniprot_entries gwas_catalog').split()

    for vcf_field in vcf_fields:
        assert vcf_field in pipeline.rs_variants
        assert vcf_field in pipeline.other_variants
    for field in annotation_fields:
        assert field in pipeline.rs_variants

    variant = pipeline.rs_variants.iloc[0]

    omim_entries = variant['omim_entries']
    assert omim_entries[0]['gene_omim_id'] == '300644'
    assert omim_entries[0]['sub_id'] == '0026'

    clinvar_entries = variant['clinvar_entries']
    assert len(clinvar_entries) == 7
    seen_rsids = {entry['rsid'] for entry in clinvar_entries}
    assert seen_rsids == {rsid}

    assert len(variant['snpeff_myvariant']) == 5

    expected_maf_keys = 'esp6500_ea_af_A 1000gp3_af_A exac_eas_af_A'.split()
    assert all(key in variant['maf'][0] for key in expected_maf_keys)

    assert variant['hgvs'][0] == {
        'genomic_allele': 'A',
        'clinvar_hgvs_c': ['LRG_672t1:c.937G>T', 'NM_000169.2:c.937G>T'],
        'clinvar_hgvs_g': ['LRG_672:g.14532G>T',
                           'NC_000023.10:g.100653420C>A',
                           'NC_000023.11:g.101398432C>A',
                           'NG_007119.1:g.14532G>T'],
        'myvariant_hgvs_g': 'chrX:g.100653420C>A',
        'snpeff_hgvs_c': ['NM_000169.2:c.937G>T',
                          'NM_021029.5:c.*475C>A',
                          'NM_001199972.1:c.*367C>A',
                          'NM_001199973.1:c.408+2975C>A',
                          'NM_001199974.1:c.285+6610C>A'],
        'snpeff_hgvs_p': ['NM_000169.2:p.Asp313Tyr']
    }

    assert variant['dbsnp_myvariant'][0]['rsid'] == rsid

    entrez_snp = variant['dbsnp_entrez']
    assert entrez_snp['rsid'] == rsid
    assert len(entrez_snp['hgvs']) == 16
    assert entrez_snp['alleles'] == 'G/T'

    uniprot_entry = variant['uniprot_entries'][0]
    assert uniprot_entry['rsid'] == rsid
    assert uniprot_entry['prot_change'] == 'p.Asp313Tyr'
    assert uniprot_entry['gene_symbol'] == 'GLA'
    assert uniprot_entry['review'] == 'in FD'

    gene_ids = {'2717', '100529097'}
    assert set(variant['entrez_gene_ids']) == gene_ids
    assert set(variant['entrez_gene_symbols']) == {'RPL36A-HNRNPH2', 'GLA'}
    assert set(pipeline.gene_annotations['entrez_gene_id'].values) == gene_ids

