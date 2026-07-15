from pathlib import Path
# ============================================================
# Project Root
# ============================================================

PROJECT_ROOT = Path(
    "/Users/anishsinha/Desktop/thesis/preprocessing"
    # "/media/rokny/DATA3/Anish/preprocessing"
)

# ============================================================
# Common Directories
# ============================================================

DATASETS_TO_PREPROCESS = (
    PROJECT_ROOT / "datasets_to_preprocess"
)

PREPROCESSED_DATA = (
    PROJECT_ROOT / "preprocessed_data"
)

GENE_ALIGNMENT = (
    PROJECT_ROOT / "gene_alignment"
)

# ============================================================
# Gene Reference Files
# ============================================================

GFF3_PATH = (
    GENE_ALIGNMENT
    / "Homo_sapiens.GRCh38.115.gff3"
)

GTF_FILE = (
    GENE_ALIGNMENT
    / "gencode.v19.chr_patch_hapl_scaff.annotation.gtf"
)