import re
import logging

from Bio import SwissProt as SwissProtReader
from Bio import ExPASy
from Bio.SeqUtils import seq3

from anotamela.annotators.base_classes import ParallelWebAnnotator


logger = logging.getLogger(__name__)


class UniprotAnnotator(ParallelWebAnnotator):
    SOURCE_NAME = 'uniprot'
    ANNOTATIONS_ARE_JSON = True

    @staticmethod
    def _query(id_):
        try:
            handle = ExPASy.get_sprot_raw(id_)
        except ValueError:
            logger.warning('"%s" not found in SwissProt' % id_)
        else:
            record = SwissProtReader.read(handle)
            # Make the record serializable for the database
            new_references = []
            for old_reference in record.references:
                new_reference = old_reference.__dict__
                new_reference['references'] = dict(old_reference.references)
                new_references.append(new_reference)
            record.references = new_references
            return record.__dict__

    @classmethod
    def _parse_annotation(cls, record):
        """
        Returns a list of variants (each represented as a dictionary) for the
        given record. All variants belong to the same protein.
        """
        prot_id = record['accessions'][0]
        prot_url = 'http://www.uniprot.org/uniprot/%s' % prot_id
        gene_match = re.search(r'Name=(.+?)( \{.+\})?;( Synonyms)?',
                               record['gene_name'])
        gene_symbol = (gene_match and gene_match.group(1)) or None

        variants = [cls._parse_variant_feature(feature)
                    for feature in record['features'] if feature[0] == 'VARIANT']
        for variant in variants:
            variant.update({
                    'prot_id': prot_id,
                    'prot_url': prot_url,
                    'gene_symbol': gene_symbol
                })
        return variants

    @staticmethod
    def _parse_variant_feature(feature):
        variant = {
            'pos': feature[1],
            'pos_stop': feature[2],
            'desc': feature[3],
            'variant_id': feature[4],
        }

        variant_url = 'http://web.expasy.org/variant_pages/{}.html'
        variant['variant_url'] = variant_url.format(variant['variant_id'])

        regex = r'(?P<old_aa>[A-Z]+) -> (?P<new_aa>[A-Z]+)'
        match = re.search(regex, variant['desc'])
        if match:
            variant['old_aa'], variant['new_aa'] = match.groups()
            variant['prot_change'] = 'p.{}{}{}'.format(seq3(variant['old_aa']),
                                                       variant['pos'],
                                                       seq3(variant['new_aa']))

        variant['pmids'] = re.findall(r'PubMed:(\d+)', variant['desc'])

        matches = re.search(r'dbSNP:(rs\d+)', variant['desc'])
        if matches:
            assert len(matches.groups()) == 1
            variant['rsid'] = matches.group(1)

        match = re.search(r'\((.+)\)', variant['desc'])
        if match:
            review = re.sub(r'; dbSNP:rs\d+', '', match.group(1))
            variant['review'] = review

        return variant

