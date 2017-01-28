from .dbsnp_web_annotator import DbsnpWebAnnotator
from .dbsnp_entrez_annotator import DbsnpEntrezAnnotator
from .dbsnp_myvariant_annotator import DbsnpMyvariantAnnotator
from .dbnsfp_annotator import DbnsfpAnnotator
from .clinvar_rs_annotator import ClinvarRsAnnotator
from .hgvs_annotator import HgvsAnnotator
from .snpeff_annotator import SnpeffAnnotator
from .maf_annotator import MafAnnotator
from .gwas_catalog_annotator import GwasCatalogAnnotator
from .ensembl_annotator import EnsemblAnnotator

from .omim_gene_annotator import OmimGeneAnnotator
from .omim_variant_annotator import OmimVariantAnnotator
from .mygene_annotator import MygeneAnnotator
from .gene_entrez_annotator import GeneEntrezAnnotator
from .uniprot_annotator import UniprotAnnotator
from .pubmed_annotator import PubmedAnnotator


RS_ANNOTATOR_CLASSES = {
    'clinvar': ClinvarRsAnnotator,
    'snpeff': SnpeffAnnotator,
    'hgvs': HgvsAnnotator,
    'maf': MafAnnotator,
    'dbsnp_myvariant': DbsnpMyvariantAnnotator,
    'dbsnp_entrez': DbsnpEntrezAnnotator,
    'dbsnp_web': DbsnpWebAnnotator,
    'dbnsfp': DbnsfpAnnotator,
    'gwas_catalog': GwasCatalogAnnotator,
    'ensembl': EnsemblAnnotator,
}

ENTREZ_GENE_ANNOTATOR_CLASSES = {
    'mygene': MygeneAnnotator,
    'entrez': GeneEntrezAnnotator,
}

