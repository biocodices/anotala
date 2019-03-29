from os.path import expanduser, abspath

from anotala.annotators.base_classes import Annotator


class LocalFileAnnotator(Annotator):
    """
    Abstract Base Class for annotators that use a local file as source of data.
    To implement a class that inherits from this one, define the methods:

        - _read_file(self, path) # => returns data
        - _parse_data(self, data) # => returns parsed data
        - _annotate_many_ids(self, ids_to_annotate) # => returns a dict with
          the annotated ids as keys and annotations as values
        - _filter_data_by(data, field_name, field_values_to_keep) [optional]

    The base class has a _load_data() method that will make use of the defined
    methods above. The loaded data is stored in a property "data".

    Then you initialize the class with the parameter *path_to_annotations_file*.
    """
    PATH_TO_ANNOTATIONS_FILE = None

    def __init__(self, path_to_annotations_file=None):
        super().__init__()
        path = path_to_annotations_file or self.PATH_TO_ANNOTATIONS_FILE
        self.path = abspath(expanduser(path))

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
