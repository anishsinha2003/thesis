import os
import pandas as pd
import scanpy as sc
import anndata as ad


def load_gse145281(dataset_dir):

    adatas = []

    for file in sorted(os.listdir(dataset_dir)):

        if not file.endswith("_raw.txt.gz"):
            continue

        sample_id = file.replace("_raw.txt.gz", "")

        path = os.path.join(dataset_dir, file)

        print(f"Loading {sample_id}")

        df = pd.read_csv(
            path,
            sep="\t",
            compression="gzip",
            index_col=0
        )

        # gene x cell -> cell x gene
        adata = sc.AnnData(df.T)

        adata.var_names_make_unique()
        adata.obs_names_make_unique()

        adata.obs["sample_id"] = sample_id

        adatas.append(adata)

    combined = ad.concat(
        adatas,
        join="outer",
        merge="same",
        fill_value=0
    )

    combined.var_names_make_unique()
    combined.obs_names_make_unique()

    return combined