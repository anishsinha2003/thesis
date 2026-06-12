# ============================================================
# Import Required Libraries
# ============================================================

# scanpy:
# Main framework used for single-cell RNA-seq analysis.
# Provides tools for preprocessing, quality control,
# dimensionality reduction, graph construction, and visualisation.

import scanpy as sc

# pandas:
# Used for reading and handling tabular files such as
# gene lists and barcode metadata.

import pandas as pd

# mmread:
# Used to load Matrix Market (.mtx) sparse matrix files,
# which store the raw gene expression counts.

from scipy.io import mmread

def load_10x_sample(sample_path):

    # ============================================================
    # Load Raw Gene Expression Matrix
    # ============================================================

    # The matrix file contains raw gene expression counts.
    # Rows originally represent genes and columns represent cells.
    # The matrix is transposed (.T) so that:
    #
    # rows   = cells
    # columns = genes
    #
    # which is the format expected by Scanpy.


    matrix = mmread(
        f"{sample_path}/matrix.mtx.gz"
    ).T

    # Convert matrix to CSR sparse format.
    # This format is more efficient for filtering,
    # indexing, and preprocessing operations.

    matrix = matrix.tocsr()

    # ============================================================
    # Load Gene Names
    # ============================================================

    # The features.tsv file contains the gene identifiers
    # corresponding to each column in the expression matrix.

    genes = pd.read_csv(
        f"{sample_path}/features.tsv.gz",
        header=None
    )

    # ============================================================
    # Load Cell Barcodes
    # ============================================================

    # The barcodes.tsv file contains unique identifiers
    # for each individual cell in the dataset.

    barcodes = pd.read_csv(
        f"{sample_path}/barcodes.tsv.gz",
        header=None
    )

    # ============================================================
    # Create AnnData Object
    # ============================================================

    # AnnData is the standard data structure used in Scanpy.
    # It stores:
    #
    # - gene expression matrix
    # - cell metadata
    # - gene metadata
    # - embeddings
    # - graphs
    # - preprocessing outputs
    #
    # in a single unified object.

    adata = sc.AnnData(matrix)

    # ============================================================
    # Assign Gene and Cell Names
    # ============================================================

    # Assign gene names to columns (variables/features).

    adata.var_names = genes[0].values

    # Assign cell barcodes to rows (observations/cells).

    adata.obs_names = barcodes[0].values

    # Print dataset dimensions:
    # n_obs  = number of cells
    # n_vars = number of genes

    return adata