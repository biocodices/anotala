from anotamela.annotators.base_classes import LocalFileAnnotator
from anotamela.helpers import ClinvarVCFParser, path_to_source_file


class ClinvarVCFAnnotator(LocalFileAnnotator):
    """
    Base class for annotators that use a local copy of the ClinVar VCF file.
    Specific Clinvar VCF Annotators will filter the entries according to
    different data, for instance, ClinvarRsVCFAnnotator gets clinvar reports
    associated with a list of RS IDs.

    To create a new ClinVar VCF Annotator that inherits from this class,
    follow the model of ClinvarRsVCFAnnotator, define:

    - SOURCE_NAME: a class variable that uniquely identifies this annotator.
    - _annotate_many_ids(self, ids_to_annotate): an instance method that loads
      the VCF data filtering by a specific field with the provided
      *ids_to_annotate*.
    """
    PACKAGED_CLINVAR_VCF = path_to_source_file('clinvar_20180729.vcf.gz')
    # Ideally, the user would provide her own clinvar vcf file, but we include
    # one in the package so that this annotator works out of the box.

    def __init__(self, path_to_annotations_file=None):
        self._parser = ClinvarVCFParser()

        if not path_to_annotations_file:
            path_to_annotations_file = self.PACKAGED_CLINVAR_VCF

        super().__init__(path_to_annotations_file)

    def __repr__(self):
        klass = self.__class__.__name__
        args = f"path_to_annotations_file='{self.path}'"
        return (f"{klass}({args})")

    def _read_file(self, path):
        return self._parser.read_file(path)

    def _parse_data(self, data):
        return self._parser.parse_data(data)

    @staticmethod
    def _filter_data_by(df, field_name, field_values_to_keep):
        filtered_df = df[df[field_name].isin(field_values_to_keep)]
        return filtered_df.reset_index(drop=True)


