from collections import defaultdict
from itertools import chain

from more_itertools import unique_everseen

from anotamela.annotators.base_classes import ParallelWebAnnotator


class DbsnpWebAnnotator(ParallelWebAnnotator):
    """
    Provider of DbSNP annotations taken from their web. Responses are cached.

        > dbsnp_web_annotator = DbsnpWebAnnotator()
        > dbsnp_web_annotator.annotate('rs123 rs268'.split())
        # => { 'rs123': ... , 'rs268': ... }

    """
    SOURCE_NAME = 'dbsnp_web'
    BATCH_SIZE = 15
    SLEEP_TIME = 1
    ANNOTATIONS_ARE_JSON = True
    MANDATORY_PROXIES = True

    @staticmethod
    def _url(rs):
        path = 'https://www.ncbi.nlm.nih.gov/projects/SNP/snp_gene.cgi?rs={0}'
        return path.format(rs)

    @classmethod
    def _parse_annotation(cls, annotation):
        annotation['rs_id'] = 'rs' + annotation['snp_id']
        # Extract this to a _parse_assemblies method
        assemblies = ['GRCh37.p13', 'GRCh38.p7']

        all_gene_models = []

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
            all_gene_models += gene_models # Will be used to extract consequences

            gene_symbols = {gene_model.get('geneSymbol')
                            for gene_model in gene_models}
            gene_symbols = [symbol for symbol in gene_symbols if symbol]
            key = '{}_gene_symbols'.format(assembly_name)
            annotation[key] = gene_symbols

            gene_ids = {gene_model.get('geneId') for gene_model in gene_models}
            gene_ids = [int(id_) for id_ in gene_ids if id_]
            key = '{}_gene_entrez_ids'.format(assembly_name)
            annotation[key] = gene_ids

        annotation['consequences_per_gene'] = \
            cls._extract_consequences_per_gene_from_gene_models(all_gene_models)

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

    @staticmethod
    def _extract_consequences_per_gene_from_gene_models(gene_models):
        """Extract the functional annotation per gene from +gene_models+.
        Returns a dictiionaty."""
        functional_annotations_per_gene = defaultdict(list)

        for gene_model in gene_models:
            gene = gene_model.get('geneSymbol', 'NO-GENE')
            alleles = gene_model.get('variation_allele', [])
            for allele in alleles:
                consequence = allele.get('fxnName')
                functional_annotations_per_gene[gene].append(consequence)

        for gene in functional_annotations_per_gene.keys():
            consequences = functional_annotations_per_gene[gene]
            unique_consequences = list(unique_everseen(consequences))
            functional_annotations_per_gene[gene] = unique_consequences

        return dict(functional_annotations_per_gene)
