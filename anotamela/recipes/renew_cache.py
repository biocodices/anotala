#!/usr/bin/env python

import os
import yaml
import datetime
from subprocess import check_output
from collections import defaultdict

from numpy import random
from more_itertools import chunked

from anotamela.annotators import ALL_ANNOTATOR_CLASSES
from anotamela.cache import create_cache


def load_mysql_credentials():
    path = os.path.expanduser("~/.mysql_credentials.yml")
    with open(path) as f:
        return yaml.load(f)


def get_table_ids_not_updated_in_n_days(table_name, n):
    threshold_date = datetime.datetime.now() - datetime.timedelta(n)
    threshold_date = threshold_date.strftime('%Y-%m-%d')

    credentials = load_mysql_credentials()
    mysql_command = f'SELECT id FROM {table_name} ' + \
                    f'WHERE {table_name}.last_updated < "{threshold_date}";'
    command = ("mysql -u{user} -p{pass} -h {host} -D {db} -e '{mysql_command}'"
               .format(**credentials, mysql_command=mysql_command))
    mysql_output = check_output(command, shell=True).decode('utf-8')
    ids = mysql_output.split('\n')[1:] # Removes the table header

    print(f'{len(ids)} IDs from this MySQL command:')
    print(mysql_command)
    print()

    return ids


def renew_cache(ids, annotator_class):
    mysql_cache = create_cache('mysql')
    proxies = {'http': 'socks5://caladan.local:9050'}

    annotator = annotator_class(cache=mysql_cache, proxies=proxies)
    failed_ids = defaultdict(list)

    # This chunking of IDs is done to avoid a single exception to ruin the
    # annotation of all IDs. We catch it instead and try to annotate the
    # chunk one ID at a time:
    chunk_size = 1_000
    for group_of_ids in chunked(ids, chunk_size):
        try:
            # This call is enough renew the cached ids:
            annotator.annotate(group_of_ids, use_cache=False, parse=False)
        except Exception:
            for id_ in group_of_ids:
                try:
                    annotator.annotate_one(id_, use_cache=False, parse=False)
                except Exception:
                    failed_ids[annotator.SOURCE_NAME].append(id_)
                    pass

    return dict(failed_ids)


def main():
    for table_name, annotator_class in ALL_ANNOTATOR_CLASSES.items():
        print('-' * 80)
        ids = get_table_ids_not_updated_in_n_days(table_name, n=30)
        # Don't reannotate all old cached entries at once, just a random subselection.
        # The idea is that the complete renewal of the cache will be done in several
        # runs of this script, not at once. The "//7" should mean "do it in seven runs",
        # and more specifically "in seven days" if this is a daily cronjob.
        subselection_of_ids = random.choice(ids, size=max(len(ids)//7, 1000))
        failed_ids = renew_cache(subselection_of_ids, annotator_class)
        print(f'{len(failed_ids)} failed ids for {table_name}: {failed_ids}')

    print('Done! Bye.')

main()
