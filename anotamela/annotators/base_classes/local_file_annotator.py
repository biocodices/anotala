from os.path import expanduser, abspath

from anotamela.annotators.base_classes import Annotator


class LocalFileAnnotator(Annotator):
    def __init__(self, path_to_annotations_file):
        super().__init__()
        self.path = abspath(expanduser(path_to_annotations_file))
        self.data = self._read_file(self.path)
        self.data = self._parse_data(self.data)

    def _read_file(self):
        raise NotImplementedError()

    def _parse_data(self):
        raise NotImplementedError()

    def _annotate_many_ids(self):
        raise NotImplementedError()
