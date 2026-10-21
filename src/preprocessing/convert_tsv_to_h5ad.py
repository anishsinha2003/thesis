
# ============================================================
# Import Required Libraries
# ============================================================

# pandas:
# Used to read and handle the tab-separated count matrix.

import pandas as pd

# anndata:
# Used to create and save AnnData objects in .h5ad format.

import anndata as ad

# scipy:
# Used to convert the count matrix into a sparse matrix,
# reducing memory usage for single-cell RNA-seq data.

from scipy.sparse import csr_matrix

from pathlib import Path


# ============================================================
# Load Single TSV Sample
# ============================================================

def load_single_tsv_sample(sample_path):

    sample_path = Path(sample_path)

    # ========================================================
    # Load Raw UMI Count Matrix
    # ========================================================

    # The TSV contains raw gene expression counts.
    # Rows originally represent genes and columns represent cells.
    #
    # The first column contains gene identifiers.
    #
    # The matrix is transposed so that:
    # rows   = cells
    # columns = genes
    #
    # This is the format expected by AnnData and Scanpy.

    df = pd.read_csv(
        sample_path,
        sep="\t",
        index_col=0
    )
    # Convert expression values to numeric
    df = df.apply(pd.to_numeric, errors="raise")

    print(f"Original matrix dimensions: {df.shape[0]} genes x {df.shape[1]} cells")

    # Transpose the matrix so cells are rows and genes are columns
    counts = df.T

    # Convert to a sparse matrix to reduce memory usage
    X = csr_matrix(counts.to_numpy())

    # Create AnnData object
    adata = ad.AnnData(X=X)

    # Set cell barcodes as observation names
    adata.obs_names = counts.index.astype(str)

    # Set gene identifiers as variable names
    adata.var_names = counts.columns.astype(str)

    # Store the raw UMI counts in a separate layer
    adata.layers["counts"] = adata.X.copy()

    # ========================================================
    # Add Basic Metadata
    # ========================================================

    adata.obs["sample_path"] = str(sample_path)
    adata.obs["sample_id"] = sample_path.stem

    # Print dataset dimensions
    print(
        f"Dataset dimensions: "
        f"{adata.n_obs} cells x {adata.n_vars} genes"
    )

    return adata


# ============================================================
# Find TSV Samples
# ============================================================

def find_tsv_samples(input_dir):
    """
    Recursively find TSV count matrix files in a directory.
    """

    input_dir = Path(input_dir)

    return sorted(
        path for path in input_dir.rglob("*.tsv")
        if path.is_file()
    )


# ============================================================
# Read Multiple TSV Samples
# ============================================================

def read_multiple_tsv_samples(input_dir):
    """
    Read multiple TSV count matrices into AnnData objects,
    then concatenate them into a single AnnData object.
    """

    adatas = []

    for sample_path in find_tsv_samples(input_dir):

        print(f"\nLoading sample: {sample_path.name}")

        adata = load_single_tsv_sample(sample_path)
        adatas.append(adata)

    if not adatas:
        raise FileNotFoundError(
            f"No TSV files found in {input_dir}"
        )

    # Concatenate all samples.
    # Genes are aligned across samples, with missing genes filled
    # with zero. Cell names are prefixed to avoid duplicates.

    adata = ad.concat(
        adatas,
        axis=0,
        join="outer",
        fill_value=0,
        label="sample_id_concat",
        keys=[a.obs["sample_id"].iloc[0] for a in adatas],
        index_unique="-"
    )

    print(
        f"\nCombined dataset dimensions: "
        f"{adata.n_obs} cells x {adata.n_vars} genes"
    )

    return adata


# ============================================================
# Save AnnData Object
# ============================================================

def save_h5ad_file(adata, output_file):
    """
    Save an AnnData object as an H5AD file.
    """

    output_file = Path(output_file)

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    adata.write_h5ad(output_file, compression="gzip")

    print(f"Saved dataset to: {output_file}")


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":

    input_file = Path(
        "/path/to/GSE166181_Raw_UMI_CountMatrix_unnormalized_.tsv"
    )

    output_file = Path(
        "/path/to/GSE166181_raw_counts.h5ad"
    )

    # Load a single TSV count matrix
    adata = load_single_tsv_sample(input_file)

    # Save the AnnData object
    save_h5ad_file(adata, output_file)
