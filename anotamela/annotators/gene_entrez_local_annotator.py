from os.path import join, dirname

import pandas as pd

from anotamela.annotators.base_classes import LocalFileAnnotator


class GeneEntrezLocalAnnotator(LocalFileAnnotator):
    """
    Annotates Entrez gene IDs or symbols using a local file downloaded from
    NCBI: ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/
    Homo_sapiens.gene_info.gz
    """
    SOURCE_NAME = 'gene_entrez_file'
    PATH_TO_ANNOTATIONS_FILE = join(dirname(dirname(__file__)),
                                    'sources',
                                    'Homo_sapiens.gene_info.2019-01-09.tsv')

    def _read_file(self, path):
        return pd.read_table(path)

    def _parse_data(self, data):
        return data # No parsing needed

    def _annotate_many_ids(self, ids_to_annotate):
        """
        Expects a list of either gene IDs or symbols (can be a mix of both).

        Returns a dataframe with data for those genes.
        """
        gene_ids = []
        gene_symbols = []

        for id_ in ids_to_annotate:
            try:
                gene_ids.append(int(id_))
            except ValueError:
                gene_symbols.append(id_)

        matching_id = self.data['GeneID'].isin(gene_ids)
        matching_symbol = self.data['Symbol'].isin(gene_symbols)

        return self.data[matching_id | matching_symbol].reset_index(drop=True)
