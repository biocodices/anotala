import logging

from concurrent.futures import (
    ProcessPoolExecutor,
    as_completed
)

from tqdm import tqdm


logger = logging.getLogger(__name__)


class Annotator():
    def __init__(self):
        self.name = self.__class__.__name__

    def annotate_one(self, id_, **kwargs):
        """Annotate one ID. Returns the annotation. Wrapper of annotate()."""
        id_ = str(id_)
        return self.annotate([id_], **kwargs).get(id_)

    def annotate(self, ids, parse=True, **kwargs):
        """
        Annotate many *ids*. Returns a dictionary where those ids are keys.

        Annotators can implement extra parsing of the data with
        _parse_annotation(annotation). This parsing is enabled by default,
        but you can disable it with parse=False and get the raw responses.
        This parsing is done for each annotation independently.

        Annotators can also implement extra parsing of the final annotations
        dictionary by defining _post_parse_annotations(annotations).
        This method is expected to receive the whole annotations dictionary.
        If parse=False, it won't be called.
        """
        ids_to_annotate = self._set_of_string_ids(ids)
        msg = '{} annotating {} ids'.format(self.name, len(ids_to_annotate))
        logger.info(msg)

        annotations = self._annotate_many_ids(ids_to_annotate, **kwargs)

        if ids_to_annotate:
            logger.info('{} found info for {}/{} ({:.2%}) IDs'
                        .format(self.name, len(annotations), len(ids_to_annotate),
                                len(annotations)/len(ids_to_annotate)))

        if parse and hasattr(self, '_parse_annotation'):
            logger.info('{} parsing {} annotations'.format(self.name,
                                                           len(annotations)))
            annotations = self._parse_annotations(annotations)
            # Sometimes, a non-empty raw response has actually no data about
            # the variant, so it generates an empty/None parsed annotation.
            # We remove those here:
            annotations = {k: v for k, v in annotations.items() if v}

        if parse and hasattr(self, '_post_parse_annotations'):
            logger.info('{} post-parsing {} annotations'.format(self.name,
                                                                 len(annotations)))
            annotations = self._post_parse_annotations(annotations)
            logger.info('Result: {} annotations'.format(len(annotations)))

        return annotations


    def _parse_annotations(self, annotations):
        """Parse a dict of annotations in parallel. Return a dictionary with
        the same keys and the parsed annotations."""
        parsed_annotations = {}
        with ProcessPoolExecutor() as executor:
            future_to_id = {}
            for id_, annotation in annotations.items():
                future = executor.submit(self._parse_annotation, annotation)
                future_to_id[future] = id_
            for future in tqdm(as_completed(future_to_id), total=len(future_to_id)):
                id_ = future_to_id[future]
                try:
                    parsed_annotations[id_] = future.result()
                except Exception as error:
                    msg = '{} parsing variant with id failed: {}'
                    raise type(error)(msg.format(self.name, id_))

        return parsed_annotations

    @staticmethod
    def _set_of_string_ids(ids):
        if isinstance(ids, str) or isinstance(ids, int):
            ids = [str(ids)]
        return set(str(id_) for id_ in set(ids) if id_)
