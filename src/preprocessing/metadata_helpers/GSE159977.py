
# ============================================================
# GSE159977 Metadata
# ============================================================

import pandas as pd


def add_metadata(adata, metadata_file):

    metadata = pd.read_csv(metadata_file)

    print("Metadata columns:")
    print(metadata.columns.tolist())

    # The cell-level metadata columns must be checked
    # before assigning patient and response labels.
    adata.obs["Cancer type"] = "Unknown"
    adata.obs["Tissue"] = "Unknown"
    adata.obs["Patient ID"] = "Unknown"
    adata.obs["Response"] = "Unknown"
    adata.obs["Sample ID"] = "Unknown"
    adata.obs["Treatment"] = "Unknown"

    # Extract GSM accession from sample_id.
    adata.obs["GSM"] = (
        adata.obs["sample_id"]
        .astype(str)
        .str.extract(r"(GSM\d+)", expand=False)
        .fillna("Unknown")
    )

    # Print the metadata so the correct clinical
    # columns can be mapped without guessing.
    print("\nMetadata preview:")
    print(metadata.head())

    return adata
