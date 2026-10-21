
import pandas as pd
import scanpy as sc
from pathlib import Path


def load_gse118056(input_dir):
    input_dir = Path(input_dir)
    file_path = input_dir / "GSE118056_raw.expMatrix.csv.gz"

    # Read the compressed count matrix.
    df = pd.read_csv(
        file_path,
        compression="gzip",
        index_col=0,
    )

    # Convert genes x cells to cells x genes.
    adata = sc.AnnData(X=df.T.to_numpy())

    adata.obs_names = df.columns.astype(str)
    adata.var_names = df.index.astype(str)

    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    return adata
