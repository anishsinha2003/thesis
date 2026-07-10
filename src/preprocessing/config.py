from pathlib import Path
from .metadata_helpers.GSE217245 import add_metadata as gse217245_metadata
from .metadata_helpers.GSE243572 import add_metadata as gse243572_metadata
from .metadata_helpers.GSE212217 import add_metadata as gse212217_metadata
from .metadata_helpers.GSE164237 import add_metadata as gse164237_metadata
from .metadata_helpers.GSE169246 import add_metadata as gse169246_metadata
from .metadata_helpers.GSE185204 import add_metadata as gse185204_metadata
from .metadata_helpers.GSE229353 import add_metadata as gse229353_metadata
# from .metadata_helpers.GSE205506 import add_metadata as gse205506_metadata
from .metadata_helpers.GSE145281 import add_metadata as gse145281_metadata



from .loaders.dataset_specfic_loaders.GSE169246 import load_gse169246
from .loaders.dataset_specfic_loaders.GSE145281 import load_gse145281
# from .metadata_helpers.GSE164237 import add_metadata as gse164237_metadata

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
    # works
    "GSE217245": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE217245",
        "output_file": PREPROCESSED_DATA / "GSE217245",
        "data_format": "10x",
        "prefix": "",
        "metadata_fn": gse217245_metadata,
    },
    # workss - memory to small
    "GSE243572": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE243572",
        "output_file": PREPROCESSED_DATA / "GSE243572",
        "data_format": "h5",
        "prefix": "",
        "metadata_fn": gse243572_metadata,
    },
    # # RDS
    # "GSE270235": {
    #     "input_dir": DATASETS_TO_PREPROCESS / "GSE270235",
    #     "output_file": PREPROCESSED_DATA / "GSE270235",
    #     "data_format": "h5",
    #     "prefix": "",
    #     "metadata_fn": gse243572_metadata,
    # },
    # works - memory to small
    "GSE212217": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE212217",
        "output_file": PREPROCESSED_DATA / "GSE212217",
        "data_format": "h5",
        "prefix": "",
        "metadata_fn": gse212217_metadata,
    },
    # cant find raw data
    # "EGAS00001004809": {
    #     "input_dir": DATASETS_TO_PREPROCESS / "EGAS00001004809",
    #     "output_file": PREPROCESSED_DATA / "EGAS00001004809",
    #     "data_format": "h5ad",
    #     "prefix": "",
    #     "metadata_fn": gse212217_metadata,
    # },
    # works - memroy too small
    "GSE169246": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE169246",
        "output_file": PREPROCESSED_DATA / "GSE169246",
        "data_format": load_gse169246,
        "prefix": "",
        "metadata_fn": gse169246_metadata,
    },
    # works
    "GSE164237": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE164237",
        "output_file": PREPROCESSED_DATA / "GSE164237",
        "data_format": "h5",
        "prefix": "",
        "metadata_fn": gse164237_metadata,
    },
    # works
    "GSE185204": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE185204",
        "output_file": PREPROCESSED_DATA / "GSE185204",
        "data_format": "10x",
        "prefix": "",
        "metadata_fn": gse185204_metadata,
    },
    # works
    "GSE229353": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE229353",
        "output_file": PREPROCESSED_DATA / "GSE229353",
        "data_format": "10x",
        "prefix": "",
        "metadata_fn": gse229353_metadata,
    },
    # # cant find responder labels
    # "GSE205506": {
    #     "input_dir": DATASETS_TO_PREPROCESS / "GSE205506",
    #     "output_file": PREPROCESSED_DATA / "GSE205506",
    #     "data_format": "10x",
    #     "prefix": "",
    #     "metadata_fn": gse205506_metadata,
    # },
    # works
    "GSE145281": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE145281",
        "output_file": PREPROCESSED_DATA / "GSE145281",
        "data_format": load_gse145281,
        "prefix": "",
        "metadata_fn": gse145281_metadata,
    },

}