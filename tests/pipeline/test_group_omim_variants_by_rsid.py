from anotamela.pipeline import group_omim_variants_by_rsid


def test_group_omim_variants_by_rsid():
    omim_variants = [
        {'rsids': ['rs1', 'rs2']},
        {'rsids': ['rs1']},
        {'rsids': ['rs2', 'rs3']},
        {}  # has no 'rsids' key, should be left out
    ]

    grouped = group_omim_variants_by_rsid(omim_variants)
    assert sorted(grouped.keys()) == ['rs1', 'rs2', 'rs3']
    assert len(grouped['rs1']) == 2
    assert len(grouped['rs2']) == 2
    assert len(grouped['rs3']) == 1

