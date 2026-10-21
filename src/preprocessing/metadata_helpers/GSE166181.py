
# ============================================================
# GSE166181 Metadata
# Metastatic Melanoma
# Immune Checkpoint Inhibitor (ICI) Therapy
# ============================================================

import pandas as pd
from src.preprocessing.config import DATASETS_TO_PREPROCESS

# ============================================================
# Metadata File
# ============================================================

METADATA_FILE = DATASETS_TO_PREPROCESS / "GSE166181_Metadata.tsv"


# ============================================================
# Helper Functions
# ============================================================

def get_patient_id(sample):

    return sample_metadata(sample, "sample")


def get_response(sample):

    return sample_metadata(sample, "response")


def get_sample_id(sample):

    return sample_metadata(sample, "orig.ident")


def get_timepoint(sample):

    return sample_metadata(sample, "subject")


def get_replicate(sample):

    return sample_metadata(sample, "replicate")


# ============================================================
# Metadata Lookup
# ============================================================

_METADATA_LOOKUP = None


def sample_metadata(cell_id, column):

    global _METADATA_LOOKUP

    if _METADATA_LOOKUP is None:
        raise RuntimeError(
            "Metadata has not been loaded. "
            "Call add_metadata(adata) first."
        )

    return _METADATA_LOOKUP.get(
        str(cell_id), {}
    ).get(column, "Unknown")


# ============================================================
# Metadata Alignment
# ============================================================

def add_metadata(adata, metadata_file=METADATA_FILE):

    global _METADATA_LOOKUP

    metadata = pd.read_csv(
        metadata_file,
        sep="\t",
        dtype={"cell_ID": str}
    )

    # Match metadata directly using cell identifiers.
    metadata = metadata.drop_duplicates(
        subset="cell_ID"
    )

    _METADATA_LOOKUP = (
        metadata.set_index("cell_ID")
        .to_dict(orient="index")
    )

    cell_ids = adata.obs_names.astype(str)

    # Check that the cell identifiers match.
    matched = cell_ids.isin(_METADATA_LOOKUP)

    print(
        f"Matched cells: {matched.sum()} / {len(cell_ids)}"
    )

    # Align metadata to the AnnData cell order.
    aligned = (
        metadata.set_index("cell_ID")
        .reindex(cell_ids)
    )

    adata.obs["Cancer type"] = (
        "Metastatic Melanoma"
    )

    adata.obs["Tissue"] = (
        "Peripheral Blood"
    )

    adata.obs["Cell type"] = (
        "CD8 T cells"
    )

    adata.obs["Patient ID"] = (
        aligned["sample"].fillna("Unknown").values
    )

    adata.obs["Sample ID"] = (
        aligned["orig.ident"].fillna("Unknown").values
    )

    adata.obs["Response"] = (
        aligned["response"].fillna("Unknown").values
    )

    adata.obs["Treatment"] = (
        "Immune Checkpoint Inhibitor (ICI)"
    )

    adata.obs["Timepoint"] = (
        aligned["subject"].fillna("Unknown").values
    )

    adata.obs["Replicate"] = (
        aligned["replicate"].fillna("Unknown").values
    )

    # Display patient-level response information.
    patient_info = (
        aligned[["sample", "response"]]
        .dropna()
        .drop_duplicates()
        .sort_values("sample")
    )

    print("\nPatient response labels:")

    for _, row in patient_info.iterrows():
        print(
            f"{row['sample']:<12} {row['response']}"
        )

    return adata
