# RS ANNOTATORS
from .dbsnp_web_annotator import DbsnpWebAnnotator
from .dbsnp_entrez_annotator import DbsnpEntrezAnnotator
from .dbsnp_myvariant_annotator import DbsnpMyvariantAnnotator
from .dbnsfp_annotator import DbnsfpAnnotator
from .clinvar_rs_annotator import ClinvarRsAnnotator
from .hgvs_annotator import HgvsAnnotator
from .snpeff_annotator import SnpeffAnnotator
from .frequencies_annotator import FrequenciesAnnotator
from .gwas_catalog_annotator import GwasCatalogAnnotator
from .ensembl_annotator import EnsemblAnnotator

from .clinvar_rcv_annotator import ClinvarRCVAnnotator
from .clinvar_variation_annotator import ClinvarVariationAnnotator

from .clinvar_rs_vcf_annotator import ClinvarRsVCFAnnotator
from .clinvar_pos_vcf_annotator import ClinvarPosVCFAnnotator

# GENE ANNOTATORS
from .mygene_annotator import MygeneAnnotator
from .gene_entrez_annotator import GeneEntrezAnnotator

from .omim_gene_annotator import OmimGeneAnnotator
from .omim_variant_annotator import OmimVariantAnnotator
from .uniprot_annotator import UniprotAnnotator
from .pubmed_annotator import PubmedAnnotator

from .biomart_regions_annotator import BiomartRegionsAnnotator


RS_ANNOTATOR_CLASSES = [
    ClinvarRsAnnotator,
    SnpeffAnnotator,
    HgvsAnnotator,
    FrequenciesAnnotator,
    DbsnpMyvariantAnnotator,
    DbsnpEntrezAnnotator,
    DbsnpWebAnnotator,
    DbnsfpAnnotator,
    GwasCatalogAnnotator,
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
