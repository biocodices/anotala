from os.path import expanduser, abspath

from anotamela.annotators.base_classes import Annotator


class LocalFileAnnotator(Annotator):
    """
    Abstract Base Class for annotators that use a local file as source of data.
    To implement a class that inherits from this one, define the methods:

        - _read_file(self, path)
        - _parse_data(self, path)
        - _annotate_many_ids(self, ids_to_annotate)
        - _filter_data_by(data, field_name, field_values_to_keep) [optional]

    """
    def __init__(self, path_to_annotations_file):
        super().__init__()
        self.path = abspath(expanduser(path_to_annotations_file))

    @property
    def data(self):
        if not hasattr(self, '_data'):
            self._data = self._load_data()
        return self._data

    def _load_data(self, filter_by={}):
        """
        Loads the data from the file at self.path to self.data.

        You can set +filter_by+ as a dictionary of {field_name: field_values}
        to filter the data by the field matches and reduce the parsing effort.
        """
        data = self._read_file(self.path)
        for field_name, field_values_to_keep in filter_by.items():
            data = self._filter_data_by(data, field_name, field_values_to_keep)
        data = self._parse_data(data)
        return data

    def _read_file(self):
        raise NotImplementedError()

    def _parse_data(self):
        raise NotImplementedError()

    def _annotate_many_ids(self):
        raise NotImplementedError()
