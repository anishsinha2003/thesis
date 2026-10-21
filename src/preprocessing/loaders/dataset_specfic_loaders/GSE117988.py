from pathlib import Path

import pandas as pd
import scanpy as sc

def load_gse117988(input_dir):
    input_dir = Path(input_dir)

    files = {
        "Tumor": "GSE117988_raw.expMatrix_Tumor.csv.gz",
        "PBMC": "GSE117988_raw.expMatrix_PBMC.csv.gz",
    }

    adatas = []

    for tissue, filename in files.items():
        file_path = input_dir / filename

        df = pd.read_csv(
            file_path,
            compression="gzip",
            index_col=0,
        )

        # Genes x cells -> cells x genes
        adata = sc.AnnData(X=df.T.to_numpy())
        adata.obs_names = df.columns.astype(str)
        adata.var_names = df.index.astype(str)

        adata.obs_names_make_unique()
        adata.var_names_make_unique()

        adata.obs["Tissue"] = tissue
        adata.obs["Patient_ID"] = "2586-4"
        adata.obs["Cancer_Type"] = "Merkel Cell Carcinoma"

        adatas.append(adata)

    # Keep genes shared by both matrices
    combined = sc.concat(
        adatas,
        join="inner",
        label="Source_Matrix",
        keys=["Tumor", "PBMC"],
        index_unique="-",
    )

    return combined
