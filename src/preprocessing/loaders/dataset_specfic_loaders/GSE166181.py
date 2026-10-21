
import os
import pandas as pd
import scanpy as sc
import anndata as ad


def load_gse166181(dataset_dir):

    file = "GSE166181_Raw_UMI_CountMatrix_unnormalized_.tsv"
    path = os.path.join(dataset_dir, file)

    print(f"Loading {file}")

    df = pd.read_csv(
        path,
        sep="\t",
        index_col=0
    )

    # gene x cell -> cell x gene
    adata = sc.AnnData(df.T)

    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    # Preserve original cell identifiers
    adata.obs["cell_id"] = adata.obs_names.astype(str)

    # Extract sample IDs from cell identifiers
    adata.obs["sample_id"] = (
        adata.obs["cell_id"]
        .str.extract(r"((?:T\d+-)?Mela\d+\s*rep\d+)", expand=False)
        .fillna("unknown")
        .str.strip()
    )

    return adata
