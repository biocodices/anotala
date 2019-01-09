from anotamela.annotators.base_classes import EntrezAnnotator
from anotamela.helpers import camel_to_snake


class GeneEntrezAnnotator(EntrezAnnotator):
    """
    Provider of Entrez Gene annotations. Usage:

        > gene_entrez = GeneEntrezAnnotator()
        > gene_entrez.annotate_one('1756')
    """
    SOURCE_NAME = 'gene_entrez'
    ANNOTATIONS_ARE_JSON = True
    ENTREZ_SERVICE = 'esummary'
    ENTREZ_PARAMS = {'db': 'gene', 'retmode': 'xml'}
    USE_ENTREZ_READER = True

    @staticmethod
    def _annotations_by_id(ids, multigene_response):
        genes = multigene_response['DocumentSummarySet']['DocumentSummary']
        id_to_gene = dict(zip(ids, genes))
        for id_, gene in id_to_gene.items():
            gene['id'] = id_
        return id_to_gene.items()

    @staticmethod
    def _parse_annotation(raw_annotation):
        fields_to_keep = ('Chromosome Description MapLocation Name Summary '
                          'NomenclatureName NomenclatureSymbol '
                          'OtherDesignations OtherAliases Mim').split()

        ann = {camel_to_snake(field): raw_annotation[field]
               for field in fields_to_keep if field in fields_to_keep}

        for key, value in raw_annotation['Organism'].items():
            new_key = 'organism_' + camel_to_snake(key).replace('i_d', 'id')
            ann[new_key] = value

        ann['url'] = 'https://www.ncbi.nlm.nih.gov/gene/' + raw_annotation['id']
        ann['entrez_id'] = int(raw_annotation['id'])
        return ann

