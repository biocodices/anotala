from anotamela.pipeline import annotate_rsids_with_clinvar


def test_annotate_rsids_with_clinvar(dict_cache):
    ids_to_annotate = [
        'rs268',  # SNP
        'rs1799990',  # Haplotype
    ]
    result = annotate_rsids_with_clinvar(ids_to_annotate,
                                         cache=dict_cache,
                                         grouped_by_rsid=True)

    # SNP:
    assert result['rs268'][0]['dbsnp_id'] == 'rs268'

    # Haplotype:
    #  assert result['rs1799990'][0]['dbsnp_ids'] == ['rs1799990', 'rs74315403']

    # Grouping by rs
    result = annotate_rsids_with_clinvar(['rs268', 'rsNonExistent'],
                                         cache=dict_cache,
                                         grouped_by_rsid=True)
    assert result['rs268'][0]['dbsnp_id'] == 'rs268'
    assert result['rsNonExistent'] == []
