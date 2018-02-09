from itertools import chain

from anotamela.annotators.base_classes import ParallelAnnotator


class DbsnpWebAnnotator(ParallelAnnotator):
    """
    Provider of DbSNP annotations taken from their web. Responses are cached.

        > dbsnp_web_annotator = DbsnpWebAnnotator()
        > dbsnp_web_annotator.annotate('rs123 rs268'.split())
        # => { 'rs123': ... , 'rs268': ... }

    """
    SOURCE_NAME = 'dbsnp_web'
    BATCH_SIZE = 15
    SLEEP_TIME = 10
    ANNOTATIONS_ARE_JSON = True

    @staticmethod
    def _url(rs):
        path = 'https://www.ncbi.nlm.nih.gov/projects/SNP/snp_gene.cgi?rs={0}'
        return path.format(rs)

    @staticmethod
    def _parse_annotation(annotation):
        annotation['rs_id'] = 'rs' + annotation['snp_id']

        # Extract this to a _parse_assemblies method
        assemblies = ['GRCh37.p13', 'GRCh38.p7']

        for assembly_name in assemblies:
            entries = annotation['assembly'].get(assembly_name) or []
            entries = [e for e in entries
                       if e['groupTerm'] == 'Primary_Assembly' and
                       e['contigLabel'] != 'PAR']
            # The 'PAR' contigLabel indicates a pseudo-autosomal region.
            # I've seen it for SNP rs6603251, located in the X chromosome
            # and in the "Y (PAR)". I decide to keep the X position.

            if not entries:
                continue

            entry = entries.pop()
            assert not entries

            key = '{}_reverse'.format(assembly_name)
            # Convert values '0' and '1' into False and True:
            annotation[key] = bool(int(entry['snp2chrOrien']))

            # Position
            keys = {
                'chr': 'chrom',
                'chrPosFrom': 'start',
                'chrPosTo': 'stop',
            }
            for key, name in keys.items():
                new_key = '{}_{}'.format(assembly_name, name)
                value = entry.get(key)
                if value and 'start' in new_key or 'stop' in new_key:
                    value = int(value)
                annotation[new_key] = value

            # Gene(s)
            gene_models = entry.get('geneModel') or []
            gene_symbols = {gene_model.get('geneSymbol')
                            for gene_model in gene_models}
            gene_symbols = [symbol for symbol in gene_symbols if symbol]
            key = '{}_gene_symbols'.format(assembly_name)
            annotation[key] = gene_symbols

        return annotation

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

