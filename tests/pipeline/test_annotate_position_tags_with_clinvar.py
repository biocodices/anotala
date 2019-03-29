from anotala.pipeline import annotate_position_tags_with_clinvar


def test_annotate_position_tags_with_clinvar(dict_cache):
    positions_to_annotate = [
        '8:19813529', # rs268
        'pos:NonExistent',
    ]
    result = annotate_position_tags_with_clinvar(positions_to_annotate,
                                                 cache=dict_cache,
                                                 grouped_by_position=True)
    assert result['8:19813529'][0]['dbsnp_id'] == 'rs268'
    assert result['rsNonExistent'] == []
