
from pathlib import Path

import anndata as ad
import pandas as pd


def load_gse115978(input_dir):

    input_dir = Path(input_dir)

    counts_file = input_dir / "GSE115978_counts.csv.gz"
    annotations_file = (
        input_dir / "GSE115978_cell.annotations.csv.gz"
    )

    # Load gene expression counts.
    counts = pd.read_csv(
        counts_file,
        index_col=0,
        compression="gzip",
    )

    # Rows are genes and columns are cells.
    counts = counts.apply(pd.to_numeric, errors="raise")

    # Transpose so rows are cells and columns are genes.
    adata = ad.AnnData(
        X=counts.T.to_numpy(),
        obs=pd.DataFrame(index=counts.columns.astype(str)),
        var=pd.DataFrame(index=counts.index.astype(str)),
    )

    adata.obs_names.name = "cell_id"
    adata.var_names.name = "gene_id"

    # Load cell annotations.
    annotations = pd.read_csv(
        annotations_file,
        compression="gzip",
    )

    # Find the column containing cell identifiers.
    possible_cell_columns = [
        "cells",
        "cell",
        "cell_id",
        "cell.name",
        "Cell",
    ]

    cell_column = next(
        (
            column
            for column in possible_cell_columns
            if column in annotations.columns
        ),
        None,
    )

    if cell_column is None:
        raise ValueError(
            "Could not find the cell ID column in "
            f"{annotations_file}. Columns: "
            f"{annotations.columns.tolist()}"
        )

    annotations[cell_column] = (
        annotations[cell_column].astype(str)
    )

    annotations = annotations.drop_duplicates(
        subset=cell_column
    ).set_index(cell_column)

    # Match annotations to the expression matrix.
    matched = annotations.reindex(adata.obs_names)

    for column in matched.columns:
        adata.obs[column] = matched[column].to_numpy()

    adata.obs["sample_id"] = (
        adata.obs["samples"].astype(str)
        if "samples" in adata.obs.columns
        else "Unknown"
    )

    adata.obs["sample_path"] = str(counts_file)

    return adata
