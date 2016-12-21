#!/usr/bin/env python
"""
Update the pickles used to test the annotators.

Usage:
    update_test_pickles.py ANNOTATOR_CLASS...

Example:
    update_test_pickles.py OmimGeneAnnotator
    update_test_pickles.py OmimGeneAnnotator DbsnpWebAnnotator

"""

import pickle
from os.path import dirname, join, abspath, basename

from docopt import docopt

from anotamela import *


def main(args):
    for annotator_class_name in args['ANNOTATOR_CLASS']:
        annotator_class = eval(annotator_class_name)
        annotator = annotator_class(cache='_dict')

        fn = pickle_path(annotator, parsed=False)
        print('Read the raw input from "{}"'.format(basename(fn)))
        with open(fn, 'rb') as f:
            raw = pickle.load(file=f)

        print('Parse it')
        parsed = annotator._parse_annotation(raw)

        fn = pickle_path(annotator, parsed=True)
        print('Dump it to "{}"'.format(basename(fn)))
        with open(fn, 'wb') as f:
            pickle.dump(obj=parsed, file=f)


def pickle_path(annotator, parsed):
    base_dir = join(dirname(dirname(dirname(abspath(__file__)))), 'tests/files')
    subdir = 'parsed_annotations' if parsed else 'raw_annotations'
    return join(base_dir, subdir, annotator.SOURCE_NAME + '.pickle')


if __name__ == '__main__':
    args = docopt(__doc__)
    main(args)

