from anotamela.pipeline import group_swissprot_variants_by_rsid


def test_group_swissprot_variants_by_rsid():
    swissprot_variants = [
        {},  # has no 'rsid'
        {'rsid': 'rs1'},
        {'rsid': 'rs2'}
    ]

    grouped = group_swissprot_variants_by_rsid(swissprot_variants)
    assert sorted(grouped.keys()) == ['rs1', 'rs2']

