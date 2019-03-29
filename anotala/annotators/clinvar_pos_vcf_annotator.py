import pandas as pd

from anotala.annotators.base_classes import ClinvarVCFAnnotator


class ClinvarPosVCFAnnotator(ClinvarVCFAnnotator):
    SOURCE_NAME = 'clinvar_pos_vcf'

    def _annotate_many_ids(self, positions_to_annotate):
        """
        Expects a list of *positions_to_annotate* as strings like "1:10000".
        It will try to find those exact positions for SNPs, and the positions
        minus 1 for deletions, in the Clinvar VCF file.
        """
        # This step is necessary to limit the expensive parsing
        # only to the VCF entries that will be actually used:
        self._data = \
            self._load_data(filter_by={'chrom_pos': positions_to_annotate})

        grouped_by_chrom_pos = self.data.groupby('chrom_pos')

        clinvar_vcf_entries_by_chrom_pos = {}

        for chrom_pos, chrom_pos_rows in grouped_by_chrom_pos:
            records = chrom_pos_rows.to_dict('records')
            clinvar_vcf_entries_by_chrom_pos[chrom_pos] = records

        return clinvar_vcf_entries_by_chrom_pos
