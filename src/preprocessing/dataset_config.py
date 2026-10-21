

from .config import DATASETS_TO_PREPROCESS
from .config import PREPROCESSED_DATA

from .metadata_helpers.GSE217245 import add_metadata as gse217245_metadata
from .metadata_helpers.GSE243572 import add_metadata as gse243572_metadata
from .metadata_helpers.GSE212217 import add_metadata as gse212217_metadata
from .metadata_helpers.GSE164237 import add_metadata as gse164237_metadata
from .metadata_helpers.GSE169246 import add_metadata as gse169246_metadata
from .metadata_helpers.GSE185204 import add_metadata as gse185204_metadata
from .metadata_helpers.GSE229353 import add_metadata as gse229353_metadata
from .metadata_helpers.GSE205506 import add_metadata as gse205506_metadata
from .metadata_helpers.GSE145281 import add_metadata as gse145281_metadata
from .metadata_helpers.GSE130157 import add_metadata as gse130157_metadata
from .metadata_helpers.GSE270235 import add_metadata as gse270235_metadata
from .metadata_helpers.GSE166181 import add_metadata as gse166181_metadata
from .metadata_helpers.GSE141383 import add_metadata as gse141383_metadata
from .metadata_helpers.GSE159977 import add_metadata as gse159977_metadata
from .metadata_helpers.GSE117988 import add_metadata as gse117988_metadata
from .metadata_helpers.GSE118056 import add_metadata as gse118056_metadata
from .metadata_helpers.GSE179784 import add_metadata as gse179784_metadata
from .metadata_helpers.SCP1493 import add_metadata as scp1493_metadata




from .metadata_helpers.bi import add_metadata as bi_metadata


from .loaders.dataset_specfic_loaders.GSE169246 import load_gse169246
from .loaders.dataset_specfic_loaders.GSE145281 import load_gse145281
from .loaders.dataset_specfic_loaders.GSE130157 import load_gse130157
from .loaders.dataset_specfic_loaders.GSE270235 import load_gse270235
from .loaders.dataset_specfic_loaders.GSE166181 import load_gse166181
from .loaders.dataset_specfic_loaders.GSE141383 import load_gse141383
from .loaders.dataset_specfic_loaders.GSE118056 import load_gse118056
from .loaders.dataset_specfic_loaders.GSE117988 import load_gse117988
from .loaders.dataset_specfic_loaders.SCP1493 import load_scp1493


from .loaders.dataset_specfic_loaders.bi import load_bi

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
        "source_assembly": "hg38",
    },

    # works - memory too small
    "GSE243572": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE243572",
        "output_file": PREPROCESSED_DATA / "GSE243572",
        "data_format": "h5",
        "prefix": "",
        "metadata_fn": gse243572_metadata,
        "source_assembly": "hg38",
    },

    # RDS
    "GSE270235": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE270235",
        "output_file": PREPROCESSED_DATA / "GSE270235",
        "data_format": load_gse270235,
        "prefix": "",
        "metadata_fn": gse270235_metadata,
        "source_assembly": "hg38",
    },

    # works - memory too small
    "GSE212217": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE212217",
        "output_file": PREPROCESSED_DATA / "GSE212217",
        "data_format": "h5",
        "prefix": "",
        "metadata_fn": gse212217_metadata,
        "source_assembly": "hg38",
    },

    # cant find raw data
    # "EGAS00001004809": {
    #     ...
    #     "source_assembly": "hg38",
    # },

    # works - memory too small
    "GSE169246": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE169246",
        "output_file": PREPROCESSED_DATA / "GSE169246",
        "data_format": load_gse169246,
        "prefix": "",
        "metadata_fn": gse169246_metadata,
        "source_assembly": "hg38",
    },

    # works
    "GSE164237": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE164237",
        "output_file": PREPROCESSED_DATA / "GSE164237",
        "data_format": "h5",
        "prefix": "",
        "metadata_fn": gse164237_metadata,
        "source_assembly": "hg38",
    },

    # works
    "GSE185204": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE185204",
        "output_file": PREPROCESSED_DATA / "GSE185204",
        "data_format": "10x",
        "prefix": "",
        "metadata_fn": gse185204_metadata,
        "source_assembly": "hg19",
    },

    # works
    "GSE229353": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE229353",
        "output_file": PREPROCESSED_DATA / "GSE229353",
        "data_format": "10x",
        "prefix": "",
        "metadata_fn": gse229353_metadata,
        "source_assembly": "hg38",
    },

    # works
    "GSE205506": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE205506",
        "output_file": PREPROCESSED_DATA / "GSE205506",
        "data_format": "10x",
        "prefix": "",
        "metadata_fn": gse205506_metadata,
        "source_assembly": "hg38",
    },

    # works
    "GSE145281": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE145281",
        "output_file": PREPROCESSED_DATA / "GSE145281",
        "data_format": load_gse145281,
        "prefix": "",
        "metadata_fn": gse145281_metadata,
        "source_assembly": "hg38",
    },

    # works
    "GSE130157": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE130157",
        "output_file": PREPROCESSED_DATA / "GSE130157",
        "data_format": load_gse130157,
        "prefix": "",
        "metadata_fn": gse130157_metadata,
        "source_assembly": "hg19",
    },

    #new datasets
    # was rds then converted to h5ad
    "GSE179784": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE179784",
        "output_file": PREPROCESSED_DATA / "GSE179784",
        "data_format": "h5ad",
        "prefix": "",
        "metadata_fn": gse179784_metadata,
        "source_assembly": "hg19",
    },

    "GSE166181": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE166181",
        "output_file": PREPROCESSED_DATA / "GSE166181",
        "data_format": "tabular",
        "prefix": "",
        "metadata_fn": gse166181_metadata,
        "source_assembly": "hg19",
    },
    # no resp lbls
    "GSE141383": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE141383",
        "output_file": PREPROCESSED_DATA / "GSE141383",
        "data_format": load_gse141383,
        "prefix": "",
        "metadata_fn": gse141383_metadata,
        "source_assembly": "hg19",
    },
    # no resp lbls
    "GSE159977": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE159977",
        "output_file": PREPROCESSED_DATA / "GSE159977",
        "data_format": "10x",
        "prefix": "",
        "metadata_fn": gse159977_metadata,
        "source_a"
        "ssembly": "hg38",
    },
    # no resp lbls
    "GSE117988" : {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE117988",
        "output_file": PREPROCESSED_DATA / "GSE117988",
        "data_format": load_gse117988,
        "prefix": "",
        "metadata_fn": gse117988_metadata,
        "source_assembly": "hg38",
    },
    # no resp lbls
    "GSE118056": {
        "input_dir": DATASETS_TO_PREPROCESS / "GSE118056",
        "output_file": PREPROCESSED_DATA / "GSE118056",
        "data_format": load_gse118056,
        "prefix": "",
        "metadata_fn": gse118056_metadata,
        "source_assembly": "GRCh38",
    },
    "SCP1493": {
        "input_dir": DATASETS_TO_PREPROCESS / "SCP1493",
        "output_file": PREPROCESSED_DATA / "SCP1493",
        "data_format": load_scp1493,
        "prefix": "",
        "metadata_fn": scp1493_metadata,
        "source_assembly": "human",
    },

}