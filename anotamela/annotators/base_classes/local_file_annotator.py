from os.path import expanduser, abspath


class LocalFileAnnotator():
    def __init__(self, path_to_annotations_file):
        self.path = abspath(expanduser(path_to_annotations_file))

    def _annotate_many_ids(self):
        raise NotImplementedError()
