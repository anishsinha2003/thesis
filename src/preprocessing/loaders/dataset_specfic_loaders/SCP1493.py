
import pandas as pd
import scanpy as sc
from pathlib import Path
from scipy.io import mmread


def load_scp1493(input_dir):
    input_dir = Path(input_dir)

    matrix_file = input_dir / "BM_c1000_raw.mtx"
    cell_file = input_dir / "BM_c1000_preQC_cellnames.tsv"
    gene_file = input_dir / "BM_c1000_preQC_genes.tsv"

    # Load the sparse gene expression matrix.
    matrix = mmread(matrix_file).tocsr()

    # Read cell and gene identifiers.
    cells = pd.read_csv(
        cell_file, sep="\t", header=None
    ).iloc[:, 0].astype(str)

    genes = pd.read_csv(
        gene_file, sep="\t", header=None
    ).iloc[:, 0].astype(str)

    # The matrix is expected to be genes x cells.
    if matrix.shape == (len(genes), len(cells)):
        matrix = matrix.T.tocsr()
    elif matrix.shape != (len(cells), len(genes)):
        raise ValueError(
            f"Matrix shape {matrix.shape} does not match "
            f"{len(cells)} cells and {len(genes)} genes."
        )

    adata = sc.AnnData(X=matrix)
    adata.obs_names = cells
    adata.var_names = genes

    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    return adata
