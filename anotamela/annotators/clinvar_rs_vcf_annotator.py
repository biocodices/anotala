from anotamela.annotators.base_classes import ClinvarVCFAnnotator


class ClinvarRsVCFAnnotator(ClinvarVCFAnnotator):
    SOURCE_NAME = 'clinvar_rs_vcf'

    def _annotate_many_ids(self, ids_to_annotate):
        """
        Expects a list of rs *ids_to_annotate* and it will try to find them
        in the ClinVar VCF file. Returns a dictionary where the rs IDs are
        the keys.
        """
        # This step is necessary to limit the expensive parsing
        # only to the VCF entries that will be actually used:
        self._data = self._load_data(filter_by={'rs_id': ids_to_annotate})

        rows_with_rs = self.data['rs_id'].isin(ids_to_annotate)
        grouped_by_rs = self.data[rows_with_rs].groupby('rs_id')

        clinvar_vcf_entries_by_rs = {}

        for rs_id, rs_rows in grouped_by_rs:
            records = rs_rows.to_dict('records')
            for record in records:
                record = {key: value for key, value in record.items() if value}
            clinvar_vcf_entries_by_rs[rs_id] = records

        return clinvar_vcf_entries_by_rs
