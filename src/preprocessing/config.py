from pathlib import Path
from .metadata_helpers.GSE217245 import add_metadata as gse217245_metadata
from .metadata_helpers.GSE243572 import add_metadata as gse243572_metadata
from .metadata_helpers.GSE212217 import add_metadata as gse212217_metadata

# ============================================================
# Project Root
# ============================================================

PROJECT_ROOT = Path(
    "/Users/anishsinha/Desktop/thesis/preprocessing"
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

# ============================================================
# Dataset Configurations
# ============================================================

DATASETS = {

    "GSE217245": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE217245",
        "output_file": PREPROCESSED_DATA / "GSE217245",
        "data_format": "10x",
        "prefix": "",
        "metadata_fn": gse217245_metadata,
    },

    "GSE243572": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE243572",
        "output_file": PREPROCESSED_DATA / "GSE243572",
        "data_format": "h5",
        "prefix": "",
        "metadata_fn": gse243572_metadata,
    },
    "GSE212217": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE212217",
        "output_file": PREPROCESSED_DATA / "GSE212217",
        "data_format": "h5",
        "prefix": "",
        "metadata_fn": gse212217_metadata,
    },

}