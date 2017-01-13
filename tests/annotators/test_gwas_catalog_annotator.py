from anotamela import GwasCatalogAnnotator


def test_gwas_catalog_annotator_parse_annotation():
    groups = []
    raw = {'grouped': {'resourcename': {'groups': groups}}}

    assert not GwasCatalogAnnotator._parse_annotation(raw)

    association_data = [
        {'rsId': ['rs123']},
        {'rsId': ['rs234']}
    ]
    study_data = []
    disease_trait_data = []

    groups.extend([
        {
            'groupValue': 'association',
            'doclist': {'docs': association_data,}
        },
        {
            'groupValue': 'study',
            'doclist': {'docs': study_data,}
        },
        {
            'groupValue': 'diseasetrait',
            'doclist': {'docs': disease_trait_data,}
        },
    ])

    parsed = GwasCatalogAnnotator._parse_annotation(raw)
    assert parsed == association_data

