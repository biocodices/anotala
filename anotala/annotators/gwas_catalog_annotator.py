import re

from anotala.annotators.base_classes import ParallelWebAnnotator
from anotala.helpers import camel_to_snake


class GwasCatalogAnnotator(ParallelWebAnnotator):
    """
    Annotates rsIDs using GWAS Catalog API:
    https://www.ebi.ac.uk/gwas/rest/docs/api#_retrieve_a_snp
    """
    SOURCE_NAME = 'gwas_catalog'
    ANNOTATIONS_ARE_JSON = True
    MANDATORY_PROXIES = True

    KEYS_TO_CAMELCASE = [
        'rsId',
        'functionalClass',
        'lastUpdateDate',
    ]

    @classmethod
    def _url(cls, id_):
        """
        Expects a dbSNP rs ID, returns the API URL for that variant.
        """
        if not re.match(r'^rs\d+$', id_):
            raise ValueError(f'Bad id: {id_}. {cls.__name__} only accepts rsIDs)')
        return f'https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{id_}'

    @classmethod
    def _parse_annotation(cls, response):
        parsed = cls._convert_keys_to_camelcase(response)
        return parsed

    @classmethod
    def _convert_keys_to_camelcase(cls, dictionary):
        result = {}
        for key in dictionary:
            if key in cls.KEYS_TO_CAMELCASE:
                result[camel_to_snake(key)] = dictionary[key]
            else:
                result[key] = dictionary[key]
        return result
