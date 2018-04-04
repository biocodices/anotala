from os.path import dirname, join, isfile


def path_to_source_file(filename):
    package_dir = dirname(dirname(__file__))
    sources_dir = join(package_dir, 'sources')
    filepath = join(sources_dir, filename)

    if not isfile(filepath):
        raise Exception(f'"{filename}" not found in {sources_dir}')

    return filepath
