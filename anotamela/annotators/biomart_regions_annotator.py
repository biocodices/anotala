from biomart import BiomartServer

from anotamela.annotators.base_classes import ParallelWebAnnotator


class BiomartRegionsAnnotator(ParallelWebAnnotator):
    """
    This class annotates a given list of regions with the SNPs located in those
    regions. GRCh37 coordinates are expected.

    Regions should be identified with the format chrom:start:end:strand.
    Strans can be 1 (forward) or -1 (reverse).

    Usage:
        > annotator = BiomartRegionsAnnotator()
        > annotator.annotate_one('1:1000000:1000100:1')

    """
    SOURCE_NAME = 'biomart_regions'
    BATCH_SIZE = 5
    SLEEP_TIME = 0.5

    ATTRIBUTES = ['refsnp_source', 'refsnp_id',
                  'chr_name', 'chrom_start', 'chrom_end', 'chrom_strand']

    def __init__(self, verbose=False, **kwargs):
        """
        Connect the annotator to Biomart's server. Set verbose=True to get info
        of the queries as they're made.
        """
        super().__init__(**kwargs)
        self.server = BiomartServer('http://grch37.ensembl.org/biomart')
        self.server.verbose = verbose
        if self.proxies and 'http' in self.proxies:
            self.server.http_proxy = self.proxies['http']
        self.database = self.server.databases['ENSEMBL_MART_SNP']
        self.human_snps = self.database.datasets['hsapiens_snp']

    def _query(self, region):
        """
        Given a single region formatted chrom:start:end:strand, return the
        raw response from BioMart with the SNPs in that region.
        """
        query_params = {'filters': {'variation_source': 'dbSNP',
                                    'chromosomal_region': region},
                        'attributes': self.ATTRIBUTES}
        response = self.human_snps.search(query_params)
        return response.content.decode()

    @classmethod
    def _parse_annotation(cls, raw_response):
        """Make a list of records (dicts) with BioMart's raw response."""
        rows = [line.split('\t') for line in raw_response.split('\n') if line]
        records = [cls._parse_row(row) for row in rows]
        return records

    @classmethod
    def _parse_row(cls, row):
        """
        Given a row like ['dbSNP', 'rs123', '1', '100', '100', '-1'],
        return a record like {'source': 'dbSNP', 'chrom_g37': '1', ... }.
        """
        names = {
            'refsnp_source': 'source',
            'refsnp_id': 'rsid',

            # Be explicit about the reference genome used:
            'chr_name': 'chrom_g37',
            'chrom_start': 'chrom_start_g37',
            'chrom_end': 'chrom_end_g37',
            'chrom_strand': 'chrom_strand_g37',
        }

        field_names = [names[attribute] for attribute in cls.ATTRIBUTES]
        record = dict(zip(field_names, row))

        record['chrom_start_g37'] = int(record['chrom_start_g37'])
        record['chrom_end_g37'] = int(record['chrom_end_g37'])
        record['chrom_strand_g37'] = int(record['chrom_strand_g37'])

        return record

