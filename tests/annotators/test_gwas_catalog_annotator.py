from anotamela import GwasCatalogAnnotator


def test_gwas_catalog_annotator():
    groups = []
    raw = {
        'grouped': {
            'resourcename': {
                'groups': groups
            }
        }
    }

    annotator = GwasCatalogAnnotator(cache='mock_cache')
    assert not annotator._parse_annotation(raw)

    association_data = [
        {'rsId': ['rs123']},
        {'rsId': ['rs234']}
    ]

    groups.extend([
        {
            'groupValue': 'association',
            'doclist': {
                'docs': association_data,
            }
        },
        {
            'groupValue': 'study',
            'doclist': {
                'docs': [
                    {'pubmedId': '123'},
                ]
            }
        },
        {
            'groupValue': 'diseasetrait',
            'doclist': {
                'docs': [
                    {'traitName': ['Pheno-1']},
                ]
            }
        },
    ])

    parsed = annotator._parse_annotation(raw)
    assert parsed == association_data

