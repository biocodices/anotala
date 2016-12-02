from os.path import isfile, expanduser

import logging

from Bio import Entrez
from tqdm import tqdm

from anotamela.annotators import AnnotatorWithCache


logger = logging.getLogger(__name__)


class ClinvarAccAnnotator(AnnotatorWithCache):
    """
    Provideer of Clinvar Accession annotations via Entrez. XML responses are
    cached and then parsed. Usage:


        > clinvar_acc_annotator = ClinvarAccAnnotator()
        > clinvar_acc_annotator.annotate(['RCV000058596'])
        # => { 'RCV000058596': ... }
    """
    SOURCE_NAME = 'clinvar_acc'

    def
