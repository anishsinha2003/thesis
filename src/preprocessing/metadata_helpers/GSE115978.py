
import os
import pandas as pd

from src.preprocessing.config import DATASETS_TO_PREPROCESS


METADATA_FILE = (
    DATASETS_TO_PREPROCESS
    / "GSE115978"
    / "GSE115978_cell.annotations.csv.gz"
)


def add_metadata(adata):
    metadata = pd.read_csv(METADATA_FILE, low_memory=False)

    print("Annotation columns:", metadata.columns.tolist())

    # Use the cell identifier column to match annotations.
    cell_column = "cells"

    if cell_column not in metadata.columns:
        raise KeyError(
            f"Expected '{cell_column}' in annotation file. "
            f"Available columns: {metadata.columns.tolist()}"
        )

    metadata[cell_column] = metadata[cell_column].astype(str).str.strip()
    metadata = metadata.drop_duplicates(subset=cell_column)
    metadata = metadata.set_index(cell_column)

    # Match annotation rows to AnnData cells.
    aligned = metadata.reindex(adata.obs_names.astype(str))

    print(
        f"Matched cells: {aligned.index.isin(metadata.index).sum()} "
        f"/ {adata.n_obs}"
    )

    # Basic metadata.
    adata.obs["Sample ID"] = (
        aligned["samples"].fillna("Unknown").astype(str).values
        if "samples" in aligned.columns
        else "Unknown"
    )

    adata.obs["Cancer type"] = "Melanoma"

    # Cell type annotation.
    if "cell.types" in aligned.columns:
        adata.obs["Cell type"] = (
            aligned["cell.types"].fillna("Unknown").astype(str).values
        )
    else:
        adata.obs["Cell type"] = "Unknown"

    # Treatment group.
    if "treatment.group" in aligned.columns:
        adata.obs["Treatment"] = (
            aligned["treatment.group"].fillna("Unknown").astype(str).values
        )
    else:
        adata.obs["Treatment"] = "Unknown"

    # Patient ID: sample names are used if no patient mapping is available.
    if "Patient" in aligned.columns:
        adata.obs["Patient ID"] = (
            aligned["Patient"].fillna("Unknown").astype(str).values
        )
    else:
        adata.obs["Patient ID"] = adata.obs["Sample ID"].values

    # This annotation file may not contain a reliable tissue column.
    adata.obs["Tissue"] = "Tumour"

    # Only assign response where treatment group supports it.
    group = adata.obs["Treatment"].astype(str).str.strip().str.lower()

    adata.obs["Response"] = group.map({
        "post-ici (resistant)": "NR",
        "or": "R",
    }).fillna("Unknown")

    return adata
