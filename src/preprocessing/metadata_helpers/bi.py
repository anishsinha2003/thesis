# src/preprocessing/metadata_helpers/bi.py

import pandas as pd

from ..config import DATASETS_TO_PREPROCESS

METADATA_FILE = DATASETS_TO_PREPROCESS / "bi/Final_SCP_Metadata.txt"


def add_metadata(adata):

    metadata = pd.read_csv(
        METADATA_FILE,
        sep="\t",
        skiprows=[1]
    )

    metadata = metadata.set_index("NAME")

    # Join metadata using cell barcodes
    adata.obs = adata.obs.join(
        metadata,
        how="left"
    )

    def map_response(response):

        if response == "ICB_PR":
            return "R"

        elif response in [
            "ICB_PD",
            "ICB_SD",
            "ICB_NE",
        ]:
            return "NR"

        return "Unknown"

    adata.obs["Cancer type"] = (
        "Clear Cell Renal Cell Carcinoma"
    )

    adata.obs["Tissue"] = (
        adata.obs["organ"]
    )

    adata.obs["Patient ID"] = (
        adata.obs["donor_id"]
    )

    adata.obs["Response"] = (
        adata.obs["ICB_Response"]
        .apply(map_response)
    )

    adata.obs["Timepoint"] = (
        adata.obs["ICB_Exposed"]
        .map({
            "ICB": "Post-Treatment",
            "No ICB": "Baseline",
        })
        .fillna("Unknown")
    )

    adata.obs["Sample ID"] = (
        adata.obs["biosample_id"]
    )

    adata.obs["Biopsy Site"] = (
        adata.obs["organ"]
    )

    adata.obs["Treatment"] = (
        "Immune Checkpoint Blockade"
    )

    adata.obs["Treatment Class"] = (
        "Anti-PD-1 / Anti-PD-L1"
    )

    print("-" * 80)

    unique_samples = (
        adata.obs[
            [
                "Sample ID",
                "Patient ID",
                "Response",
            ]
        ]
        .drop_duplicates()
        .sort_values("Sample ID")
    )

    for _, row in unique_samples.iterrows():
        print(
            f"{row['Sample ID']:<15} "
            f"{row['Patient ID']:<8} "
            f"{row['Response']}"
        )

    print("-" * 80)

    return adata