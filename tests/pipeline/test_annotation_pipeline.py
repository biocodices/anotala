import pickle

import pytest

from anotamela.cache import Cache, DictCache
from anotamela.pipeline import AnnotationPipeline
from anotamela.annotators.base_classes.parallel_web_annotator import NoProxiesException


def test_init():
    pipe = AnnotationPipeline(cache='mock_cache', proxies={}, sleep_time=1,
                              use_cache=True, use_web=False)

    # Test it creates a Cache instance for the annotation
    assert isinstance(pipe.annotation_kwargs['cache'], Cache)

    # Test it accepts an already instantiated Cache
    dict_cache = DictCache()
    pipe = AnnotationPipeline(cache=dict_cache, proxies={})
    assert pipe.annotation_kwargs['cache'] is dict_cache

    # Test it raises if it's initialized withoug explicit proxies setting
    with pytest.raises(NoProxiesException):
        AnnotationPipeline(cache=dict_cache)


def test_run_from_rsids(proxies):
    dict_cache = DictCache()

    with open(pytest.helpers.file('dict_cache_storage.pickle'), 'rb') as f:
        dict_cache.storage = pickle.load(f)

    # Every annotation will be brought from the pickle-loaded cache
    pipe = AnnotationPipeline(cache=dict_cache, use_web=False, proxies=proxies)

    # The .pickle loaded above has annotations for this single rs ID,
    # to allow this test to run fast with real data.
    rsid = 'rs1799983'
    # The pickle file was the result of pickle.dump()ing the
    # pipe.annotation_kwargs['cache'].storage after choosing 'dict' as the
    # cache of an annotation and pipe.run_from_rsids(['rs1799983']).
    # The chosen SNP has info from ClinVar, OMIM and Uniprot, but not from
    # GWAS Catalog, so that's a blind spot in this test.

    pipe.run_from_rsids([rsid])

    assert len(pipe.rs_variants) == 1
    assert len(pipe.gene_annotations) == 1

    annotations_to_expect = [
        'id',
        'dbsnp_entrez',
        'dbsnp_myvariant',
        'frequencies',
        'hgvs',
        'clinvar_entries',
        'clinvar_vcf_entries',
        'clinvar_variations',
        'clinvar_variations_from_position',
        # 'gwas_catalog',
        'dbnsfp',
        'snpeff_myvariant',
        'entrez_gene_ids',
        'entrez_gene_symbols',
        'omim_entries',
        'uniprot_entries',
        'ensembl',
    ]
    for field in annotations_to_expect:
        assert field in pipe.rs_variants.columns

    for field in ['entrez_gene_id', 'mygene', 'gene_entrez']:
        assert field in pipe.gene_annotations

    # Some small checks to kinda guess that the data was munged correctly
    variant = pipe.rs_variants.iloc[0]

    clinvar_entry = variant['clinvar_entries'][0]
    assert clinvar_entry['variant_id'] == 14015

    gene = pipe.gene_annotations.iloc[0]
    assert gene['entrez_gene_id'] == '4846'

    uniprot_entry = variant['uniprot_entries'][0]
    assert uniprot_entry['variant_id'] == 'VAR_008037'

    omim_entry = variant['omim_entries'][0]
    assert omim_entry['variant_id'] == '163729.0001'

    entrez_snp = variant['dbsnp_entrez']
    assert entrez_snp['alleles'] == 'G/T'

    myvariant_snp = variant['dbsnp_myvariant'][0]
    assert myvariant_snp['vartype'] == 'snp'

    hgvs = variant['hgvs'][0]
    assert hgvs['myvariant_hgvs_g'] == 'chr7:g.150696111T>G'

    freqs = variant['frequencies']
    assert freqs['T']['dbSNP']['General'] == 0.1763
    assert freqs['G']['dbNSFP_ExAC']['Finnish (FIN)'] == 0.7079

    dbnsfp = variant['dbnsfp'][0]
    assert dbnsfp['cds_strand'] == '+'

    snpeff = variant['snpeff_myvariant'][0]
    assert snpeff['feature_type'] == 'transcript'

    assert variant['entrez_gene_symbols'] == ['NOS3']
    assert variant['entrez_gene_ids'] == ['4846']
