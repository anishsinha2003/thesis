from pathlib import Path

import anndata as ad
import pandas as pd

from scipy.io import mmread


def load_gse169246(input_dir):

    matrix_file = input_dir / "allSamplesCombined/matrix.mtx.gz"
    features_file = input_dir / "allSamplesCombined/features.tsv.gz"
    barcodes_file = input_dir / "allSamplesCombined/barcodes.tsv.gz"

    print("Loading matrix...")
    X = mmread(matrix_file).T.tocsr()

    print("Loading features...")
    genes = pd.read_csv(
        features_file,
        sep="\t",
        header=None,
    )

    print("Loading barcodes...")
    barcodes = pd.read_csv(
        barcodes_file,
        sep="\t",
        header=None,
    )


    adata = ad.AnnData(X)

    adata.var_names = genes[0].astype(str)
    adata.obs_names = barcodes[0].astype(str)

    adata.var_names_make_unique()

    adata.obs["sample_path"] = str(input_dir)
    adata.obs["sample_id"] = "GSE169246"

    print(
        f"Dataset dimensions: "
        f"{adata.n_obs} cells x {adata.n_vars} genes"
    )

    return adata