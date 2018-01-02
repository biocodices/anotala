from anotamela import DbsnpEntrezAnnotator
#  import pytest

#  @pytest.fixture
#  def raw_annotation():
    #  return """
#  <rs bitfield="050268000A050405133F0100" rsid="234">

    #  <mergehistory rsid="123"></mergehistory>

#  </rs>
#  """


def test_post_parse_annotations():
    parsed_annotations = {
        'rs234': {
            'rsid': 'rs234',
            'synonyms': ['rs123', 'rs345'],
        }
    }

    result = DbsnpEntrezAnnotator._post_parse_annotations(parsed_annotations)

    # Keep the new ID
    assert result['rs234'] == parsed_annotations['rs234']
    # Add the synonym IDs
    assert result['rs123'] == parsed_annotations['rs234']
    assert result['rs345'] == parsed_annotations['rs234']

