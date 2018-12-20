from anotamela.annotators import *

RS_ANNOTATOR_CLASSES = [
    ClinvarRsAnnotator,
    SnpeffAnnotator,
    HgvsAnnotator,
    FrequenciesAnnotator,
    DbsnpMyvariantAnnotator,
    # DbsnpEntrezAnnotator,
    DbsnpWebAnnotator,
    DbnsfpAnnotator,
    #  GwasCatalogAnnotator,
    #  GwasCatalogOldAnnotator,
    EnsemblAnnotator,
]

RS_ANNOTATOR_CLASSES = {klass.SOURCE_NAME: klass
                        for klass in RS_ANNOTATOR_CLASSES}

ENTREZ_GENE_ANNOTATOR_CLASSES = [
    MygeneAnnotator,
    GeneEntrezAnnotator,
]

ENTREZ_GENE_ANNOTATOR_CLASSES = {klass.SOURCE_NAME: klass
                                 for klass in ENTREZ_GENE_ANNOTATOR_CLASSES}

# This is used elsewhere to programatically update the cache for
# all annotators:
ALL_ANNOTATOR_CLASSES = [
    ClinvarRCVAnnotator,
    ClinvarRsVCFAnnotator,
    ClinvarPosVCFAnnotator,
    ClinvarVariationAnnotator,
    OmimGeneAnnotator,
    OmimVariantAnnotator,
    UniprotAnnotator,
    PubmedAnnotator,
    BiomartRegionsAnnotator,
]
ALL_ANNOTATOR_CLASSES = {klass.SOURCE_NAME: klass
                         for klass in ALL_ANNOTATOR_CLASSES}
ALL_ANNOTATOR_CLASSES.update(RS_ANNOTATOR_CLASSES)
ALL_ANNOTATOR_CLASSES.update(ENTREZ_GENE_ANNOTATOR_CLASSES)
