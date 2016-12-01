import requests
from itertools import chain

from anotamela.annotators import AnnotatorWithCache


class DbsnpWebAnnotator(AnnotatorWithCache):
    """
    Provider of DbSNP annotations taken from their web. Responses are cached.

        > dbsnp_web_annotator = DbsnpWebAnnotator()
        > dbsnp_web_annotator.annotate('rs123 rs268'.split())
        # => { 'rs123': ... , 'rs268': ... }

    """
    SOURCE_NAME = 'dbsnp_web'

    @staticmethod
    def _url(rs):
        path = 'https://www.ncbi.nlm.nih.gov/projects/SNP/snp_gene.cgi?rs={0}'
        return path.format(rs)

    def _query(self, rs):
        """
        Query NCBI's dbSNP site for a given rs ID. Returns a dict or None.
        """
        url = self._url(rs)
        headers = {'Content-Type': 'application/json', 'Accept': 'application/json'}
        response = requests.get(url, headers)

        if response.ok:
            return response.json()
        else:
            print('{} status code for "{}"'.format(response.status_code, rs))
            return None

    def genes(self, rs, use_web=False):
        """Annotate the genes for a given rs."""
        ann = self.annotate(rs, use_web=use_web).get(rs)
        if not ann:
            return []

        mappings = ann.get('assembly', {}).values()
        mappings = chain(*mappings)

        gene_models = [m['geneModel'] for m in mappings if 'geneModel' in m]
        gene_models = chain(*gene_models)

        unique_genes = set()
        for gene in gene_models:
            gene_str = gene['geneSymbol'] or ''
            genes_set = set(gene_str.split('|'))
            unique_genes.update(genes_set)

        return sorted(list(unique_genes))

