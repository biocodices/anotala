import re
import logging

import pandas as pd
from more_itertools import collapse
from vcf_to_dataframe import vcf_to_dataframe


logger = logging.getLogger(__name__)


class ClinvarVCFParser:
    def read_file(self, path):
        df = vcf_to_dataframe(path)

        # This step is done before the rest of the parsing to expose the rs_id,
        # so that the data can be filtered by those IDs before the
        # computationally heavy parsing done at `parse_data`:
        df = self._add_rs_to_dbsnp_id(df)
        df['rs_id'] = df['info'].map(lambda d: d.get('rs_id'))

        return df

    def parse_data(self, df):
        logger.info(f'Parsing ClinVar VCF, this might take a while')
        df = self._extract_keys(df)
        df = self._extract_external_db_disease_ids(df)
        df = self._parse_and_extract_sources(df)
        df = self._parse_and_extract_gene_info(df)
        df = self._parse_origin(df)
        df = self._extract_molectular_consequences(df)
        df = self._rename_columns(df)
        df = self._drop_columns(df)
        return df

    @staticmethod
    def _add_rs_to_dbsnp_id(df):
        # NOTE: This method modifies the dictionaries of the original df.
        # We don't really care here, but it's not nice.
        for info_dict in df['info'].values:
            if 'RS' in info_dict:
                info_dict['rs_id'] = 'rs' + info_dict['RS']
        return df

    @staticmethod
    def _extract_keys(df):
        """
        Extracts all keys under "info" dictionaries into new columns.
        """
        df = df.copy()
        all_keys_seen = {key for info in df['info'].values
                         for key in info.keys()}
        for key in all_keys_seen:
            df[key] = df['info'].map(lambda info: info.get(key))
        return df

    @staticmethod
    def _extract_external_db_disease_ids(df):
        df = df.copy()

        def extract_disease_ids(row):
            disease_ids = []
            db_names_and_ids = row['info'].get('CLNDISDB')

            db_names_and_ids = db_names_and_ids

            if db_names_and_ids:
                # Turn "MedGen:1|MedGen:2,OMIM:3" into "MedGen:1,MedGen2,OMIM:3"
                # to ease the parsing:
                db_names_and_ids = db_names_and_ids.replace('|', ',')
                disease_ids += db_names_and_ids.split(',')

            return disease_ids

        df['disease_ids'] = df.apply(extract_disease_ids, axis=1)

        return df

    @staticmethod
    def _parse_and_extract_sources(df):
        df = df.copy()
        regex = re.compile(r'(.+?):(.+?)(?:,|$)')

        def parse_CLNVI_field(clnvi):
            name_value_pairs = regex.findall(clnvi)
            return dict(name_value_pairs)

        df['sources'] = df['CLNVI'].map(parse_CLNVI_field, na_action='ignore')

        sources_to_extract = {
            'OMIM_Allelic_Variant': 'omim_variant_id',
            'PharmGKB': 'pharmgkb_id',
            'UniProtKB_(protein)': 'uniprot_variant_id',
        }
        for source, field_name in sources_to_extract.items():
            df[field_name] = df['sources'].map(lambda d: d.get(source),
                                               na_action='ignore')
        return df

    @staticmethod
    def _parse_and_extract_gene_info(df):
        df = df.copy()

        def parse_geneinfo(geneinfo):
            info = {}
            if geneinfo and not pd.isnull(geneinfo):
                info = dict([key_val.split(':')
                             for key_val in geneinfo.split('|')])
                info = {symbol: int(id_) for symbol, id_ in info.items()}
            return info

        gene_dicts = df['GENEINFO'].map(parse_geneinfo)
        df['genes'] = gene_dicts
        df['gene_symbols'] = gene_dicts.map(lambda d: list(d.keys()))
        df['gene_entrez_ids'] = gene_dicts.map(lambda d: list(d.values()))
        return df

    @staticmethod
    def _extract_molectular_consequences(df):
        df = df.copy()

        def parse_mc_value(mc):
            if not mc or pd.isnull(mc):
                return []
            consequences = [value.split(',') for value in mc.split('|')]
            consequences = [c for c in collapse(consequences) if not 'SO:' in c]
            return consequences

        df['molecular_consequences'] = df['MC'].map(parse_mc_value)
        return df

    @staticmethod
    def _parse_origin(df):
        df = df.copy()
        origin_names = {
            '0': 'unknown',
            '1': 'germline',
            '2': 'somatic',
            '4': 'inherited',
            '8': 'paternal',
            '16': 'maternal',
            '32': 'de-novo',
            '64': 'biparental',
            '128': 'uniparental',
            '256': 'non-tested',
            '512': 'tested-inconclusive',
            '1073741824': 'other',
        }
        df['origin'] = df['ORIGIN'].map(
            lambda key: origin_names.get(key, key),
            na_action='ignore'
        )
        return df

    @staticmethod
    def _drop_columns(df):
        columns_to_drop = [
            'RS',
            'CLNVI',
            'GENEINFO',
            'ORIGIN',
            'MC',
        ]
        return df.copy().drop(columns_to_drop, axis=1)

    @staticmethod
    def _rename_columns(df):
        column_names = {
            'id': 'variation_id',
            'ALLELEID': 'allele_id',
            'CLNHGVS': 'clinvar_hgvs',
            'CLNDN': 'disease_name',
            'CLNREVSTAT': 'revision_status',
            'info': 'source_info',
            'CLNSIG': 'clinical_significance',
            'CLNVC': 'variant_type',
            'CLNVCSO': 'sequence_ontology_id',
        }
        return df.rename(columns=column_names)
