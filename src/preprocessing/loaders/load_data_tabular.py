
# Read tabular data file into an AnnData object (TSV, CSV, TXT)

import anndata as ad
import pandas as pd
import numpy as np


def guess_sep(path) -> str:
    path = str(path).lower()
    if path.endswith(".csv") or path.endswith(".csv.gz"):
        return ","
    return "\t"


def read_tabular_file(file_path):
    sep = guess_sep(file_path)

    df = pd.read_csv(
        file_path,
        sep=sep,
        index_col=0,
        low_memory=False
    )

    # Convert expression values to numeric.
    # Raise an error if any values cannot be converted.
    df = df.apply(pd.to_numeric, errors="raise")

    # Check that the expression matrix contains valid values.
    if df.isna().any().any():
        raise ValueError(
            f"Missing values found in expression matrix: {file_path}"
        )

    # Create AnnData object.
    adata = ad.AnnData(
        X=df.to_numpy(dtype=np.float32)
    )

    adata.obs_names = df.index.astype(str)
    adata.var_names = df.columns.astype(str)
    adata.var_names_make_unique()

    adata.obs["source_count_file"] = str(file_path)
    adata.obs["sample_id"] = (
        file_path.stem
        .replace(".tsv", "")
        .replace(".csv", "")
        .replace(".txt", "")
    )
    adata.obs["sample_path"] = str(file_path)

    adata.uns["input_format"] = "tabular"
    adata.uns["input_count_files"] = [str(file_path)]

    return adata


def find_count_tables(input_dir):
    """Find tabular count files in a folder."""

    patterns = [
        "*.csv", "*.tsv", "*.txt",
        "*.csv.gz", "*.tsv.gz", "*.txt.gz"
    ]

    files = []

    for pattern in patterns:
        files.extend(input_dir.rglob(pattern))

    files = sorted(files)

    if not files:
        raise FileNotFoundError(
            f"No tabular count files found in: {input_dir}"
        )

    return files


def load_tabular_folder(input_dir):
    """Load and concatenate tabular expression files."""

    tabular_files = find_count_tables(input_dir)

    adatas = [
        read_tabular_file(file_path)
        for file_path in tabular_files
    ]

    return ad.concat(
        adatas,
        axis=0,
        join="outer",
        fill_value=0,
        label="sample_id_concat",
        keys=[
            a.obs["sample_id"].iloc[0]
            for a in adatas
        ],
        index_unique=None
    )
