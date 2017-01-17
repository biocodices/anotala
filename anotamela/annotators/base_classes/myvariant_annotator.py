from operator import itemgetter
from itertools import groupby, chain
import logging

from myvariant import MyVariantInfo
from tqdm import tqdm

from anotamela.annotators.base_classes import AnnotatorWithCache
from anotamela.helpers import (
    grouped,
    infer_annotated_allele,
)


logger = logging.getLogger(__name__)


class MyVariantAnnotator(AnnotatorWithCache):
    """
    This class is meant as a class for annotators that use myvariant.info.
    To use it, define a new class that inherits from this one and that has:

        - SOURCE_NAME class variable
        - SCOPES class variable (scopes to query with the passed IDs)
        - FIELDS class variable (fields to retrieve for each ID)
        - optionally, a @staticmethod or @classmethod _parse_hit(hit) to parse
          each of the hits (results) associated to a single ID. If not defined,
          the results won't be parsed in any way.
        - optionally, set a class variable VERBOSE to get all the output
          produced by myvariant (silenced by default)
    """
    ANNOTATIONS_ARE_JSON = True
    VERBOSE = False
    BATCH_SIZE = 1000

    def _batch_query(self, ids):
        """
        Uses myvariant.info service to query many IDs across the scopes defined
        by the class variable SCOPES. It returns a dict of {id: result, ... }
        with the IDs that were found (i.e. leaves out the not found ones) and
        info from the fields defined in the class variable FIELDS.

        The successful results are cached.
        """
        if not hasattr(self, 'mv'):
            self.mv = MyVariantInfo()

        grouped_ids = list(grouped(ids, self.BATCH_SIZE))
        for batch_of_ids in tqdm(grouped_ids):
            logger.debug('{} query {} IDs'.format(self.name, len(batch_of_ids)))
            hits = self.mv.querymany(batch_of_ids, scopes=self.SCOPES,
                                     fields=self.FIELDS, verbose=self.VERBOSE)

            batch_annotations = {}
            for query, hits in groupby(hits, itemgetter('query')):
                hits = list(hits)
                if 'notfound' not in hits[0]:
                    batch_annotations[query] = hits

            yield batch_annotations

    @classmethod
    def _parse_annotation(cls, hits):
        # Gather the hits (i.e. different alleles) in the same order always:
        hits = sorted(hits, key=itemgetter('_id'))

        for hit in hits:
            hit['allele'] = infer_annotated_allele(hit['_id'])

        annotations = [cls._parse_hit(hit) for hit in hits]
        annotations = [ann for ann in annotations if ann]

        # If the annotations are already lists, merge them into a flat list:
        if all(isinstance(ann, list) for ann in annotations):
            annotations = list(chain.from_iterable(annotations))

        if annotations:
            return annotations

    @staticmethod
    def _parse_hit(hit):
        return hit

