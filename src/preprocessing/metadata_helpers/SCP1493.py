
# ============================================================
# SCP1493 Metadata
# Melanoma Brain Metastases
# Immune Checkpoint Inhibitor (ICI)
# ============================================================

import pandas as pd
from pathlib import Path
from ..config import DATASETS_TO_PREPROCESS


# ------------------------------------------------------------
# Response
# ------------------------------------------------------------

def get_response(response):

    response = str(response).strip().lower()

    if response in ["responder", "response", "responders"]:
        return "R"

    elif response in [
        "non-responder",
        "non responder",
        "nonresponders",
        "progression",
    ]:
        return "NR"

    elif response in [
        "partial responder",
        "partial response",
        "partial responders",
    ]:
        return "Partial"

    return "Unknown"


# ------------------------------------------------------------
# Metadata
# ------------------------------------------------------------

def add_metadata(adata):

    input_dir = Path(DATASETS_TO_PREPROCESS) / "SCP1493"

    metadata = pd.read_csv(
        input_dir / "BM_c1000_base_metadata.tsv",
        sep="\t",
    )

    if "NAME" not in metadata.columns:
        raise ValueError(
            "Metadata has no NAME column. "
            f"Available columns: {metadata.columns.tolist()}"
        )

    metadata["NAME"] = metadata["NAME"].astype(str)

    metadata = (
        metadata.drop_duplicates("NAME")
        .set_index("NAME")
    )

    # --------------------------------------------------------
    # Match metadata to cells
    # --------------------------------------------------------

    matched = metadata.reindex(adata.obs_names.astype(str))

    for column in matched.columns:
        adata.obs[column] = matched[column].to_numpy()

    missing = matched.isna().all(axis=1)

    if missing.any():
        print(
            f"Warning: metadata was not matched for "
            f"{missing.sum()} cells."
        )

    # --------------------------------------------------------
    # Standard metadata fields
    # --------------------------------------------------------

    adata.obs["Cancer type"] = "Melanoma"
    adata.obs["Cancer Subtype"] = "Melanoma"
    adata.obs["Disease Stage"] = "Metastatic"
    adata.obs["Tissue"] = "Brain metastasis"

    adata.obs["Treatment"] = (
        "Immune Checkpoint Inhibitor"
    )

    adata.obs["Treatment Class"] = (
        "Immune Checkpoint Inhibitor"
    )

    # --------------------------------------------------------
    # Sample ID
    # --------------------------------------------------------

    if "biosample_id" in adata.obs.columns:
        adata.obs["Sample ID"] = (
            adata.obs["biosample_id"].astype(str)
        )
    else:
        adata.obs["Sample ID"] = "Unknown"

    # --------------------------------------------------------
    # Patient ID
    # --------------------------------------------------------

    if "donor_id" in adata.obs.columns:
        adata.obs["Patient ID"] = (
            adata.obs["donor_id"].astype(str)
        )
    else:
        adata.obs["Patient ID"] = "Unknown"

    # --------------------------------------------------------
    # Response
    # --------------------------------------------------------

    # Only map response if the relevant column exists.
    # The correct source column still needs to be confirmed.

    response_columns = [
        "response",
        "Response",
        "clinical_response",
    ]

    response_column = next(
        (
            column
            for column in response_columns
            if column in adata.obs.columns
        ),
        None,
    )

    if response_column is not None:
        adata.obs["Response"] = (
            adata.obs[response_column].apply(get_response)
        )
    else:
        adata.obs["Response"] = "Unknown"

    return adata
