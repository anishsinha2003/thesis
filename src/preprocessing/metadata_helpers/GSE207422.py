
import pandas as pd

from src.preprocessing.config import DATASETS_TO_PREPROCESS


METADATA_FILE = (
    DATASETS_TO_PREPROCESS
    / "GSE207422"
    / "GSE207422_NSCLC_scRNAseq_metadata.xlsx"
)


def add_metadata(adata):
    metadata = pd.read_excel(METADATA_FILE)

    # Clean sample names in the metadata.
    metadata["Sample"] = (
        metadata["Sample"].astype(str).str.strip()
    )

    # Cell IDs look like BD_immune01_612637.
    # Extract BD_immune01 to match the Excel metadata.
    cell_ids = adata.obs_names.astype(str)

    sample_ids = cell_ids.str.extract(
        r"^(BD_immune\d+)_",
        expand=False
    )

    adata.obs["Sample ID"] = sample_ids.fillna("Unknown").values

    # Match each cell to its original sample.
    metadata_lookup = metadata.drop_duplicates(
        subset="Sample"
    ).set_index("Sample")

    matched = adata.obs["Sample ID"].map(metadata_lookup["Patient"])

    adata.obs["Patient ID"] = matched.fillna("Unknown").values

    # Add response labels from RECIST.
    recist_lookup = metadata_lookup["RECIST"]

    recist = (
        adata.obs["Sample ID"]
        .map(recist_lookup)
        .astype("string")
        .str.strip()
        .str.upper()
    )

    adata.obs["Response"] = (
        recist.map({
            "CR": "R",
            "PR": "R",
            "SD": "NR",
            "PD": "NR",
        })
        .fillna("Unknown")
        .values
    )

    # Add the remaining metadata fields.
    for column in ["Sex", "Age"]:
        adata.obs[column] = (
            adata.obs["Sample ID"]
            .map(metadata_lookup[column])
            .fillna("Unknown")
            .values
        )

    adata.obs["Cancer type"] = "NSCLC"
    adata.obs["Tissue"] = "Unknown"

    adata.obs["Treatment"] = "Toripalimab + Chemotherapy"

    print("Metadata matching:")
    print(
        "Matched samples:",
        adata.obs["Sample ID"].nunique()
    )
    print(
        "Unknown patient IDs:",
        (adata.obs["Patient ID"] == "Unknown").sum()
    )
    print("\nResponse counts:")
    print(adata.obs["Response"].value_counts())

    return adata
