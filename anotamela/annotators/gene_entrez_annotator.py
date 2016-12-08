from anotamela.annotators.base_classes import EntrezAnnotator
from anotamela.helpers import camel_to_snake


class GeneEntrezAnnotator(EntrezAnnotator):
    """
    Provider of Entrez Gene annotations. Usage:

        > gene_entrez = GeneEntrezAnnotator()
        > gene_entrez.annotate('1756')  # Also accepts a list of IDs
        # => { '1756': ... }
    """
    SOURCE_NAME = 'gene_entrez'
    ANNOTATIONS_ARE_JSON = True
    ENTREZ_PARAMS = {'db': 'gene', 'retmode': 'xml', 'service': 'esummary'}
    USE_ENTREZ_READER = True

    @staticmethod
    def _annotations_by_id(ids, multigene_response):
        genes = multigene_response['DocumentSummarySet']['DocumentSummary']
        yield from dict(zip(ids, genes)).items()

    @staticmethod
    def _parse_annotation(raw_annotation):
        fields_to_keep = ('Chromosome Description MapLocation Name Summary '
                          'NomenclatureName NomenclatureSymbol '
                          'OtherDesignations OtherAliases').split()

        ann = {camel_to_snake(field): raw_annotation[field]
               for field in fields_to_keep}

        for key, value in raw_annotation['Organism'].items():
            new_key = 'organism_' + camel_to_snake(key).replace('i_d', 'id')
            ann[new_key] = value

        return ann

