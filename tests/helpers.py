from os.path import join, dirname


def get_test_file(filename):
    return join(dirname(__file__), 'files', filename)

