import logging

from vcf_to_dataframe import vcf_to_dataframe

from anotamela.annotators.base_classes import LocalFileAnnotator


logger = logging.getLogger(__name__)


class ClinvarRsVCFAnnotator(LocalFileAnnotator):
    SOURCE_NAME = 'clinvar_rs_vcf'

    def __init__(self, path_to_annotations_file):
        super().__init__(path_to_annotations_file)
        self.data = self._read_vcf(self.path)
        self._parse_dataframe(self.data)

    @staticmethod
    def _read_vcf(path):
        return vcf_to_dataframe(path)

    @staticmethod
    def _parse_dataframe(dataframe):
        all_keys_seen = {key for info in dataframe['info'].values
                         for key in info.keys()}

        for key in all_keys_seen:
            dataframe[key] = dataframe['info'].map(lambda info: info.get(key))

        colnames = {
            'RS': 'rs_id',
            'CLNHGVS': 'clinvar_name',
        }
        dataframe.rename(columns=colnames, inplace=True)

    def _annotate_many_ids(self, ids_to_annotate):
        pass
