from anotamela.annotators.base_classes import LocalFileAnnotator
from anotamela.helpers import ClinvarVCFParser, path_to_source_file


class ClinvarRsVCFAnnotator(LocalFileAnnotator):
    SOURCE_NAME = 'clinvar_rs_vcf'

    def __init__(self, path_to_annotations_file=None):
        self._parser = ClinvarVCFParser()

        if not path_to_annotations_file:
            path_to_annotations_file = path_to_source_file('clinvar_20180128.vcf.gz')

        super().__init__(path_to_annotations_file)

    def _read_file(self, path):
        return self._parser.read_file(path)

    def _parse_data(self, data):
        return self._parser.parse_data(data)

    @staticmethod
    def _filter_data_by(df, field_name, field_values_to_keep):
        filtered_df = df[df[field_name].isin(field_values_to_keep)]
        return filtered_df.reset_index(drop=True)

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

        records_by_rs = {}

        for rs_id, rs_rows in grouped_by_rs:
            records = rs_rows.to_dict('records')
            for record in records:
                record = {key: value for key, value in record.items() if value}
            records_by_rs[rs_id] = records

        return records_by_rs
