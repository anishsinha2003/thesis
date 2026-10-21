
import os
import pandas as pd
import scanpy as sc
import anndata as ad


def load_gse141383(dataset_dir):

    adatas = []

    for file in sorted(os.listdir(dataset_dir)):

        if not file.endswith(".counts.txt.gz"):
            continue

        # Example: GSM4202145_PJ053.counts.txt.gz
        sample_id = file.replace(".counts.txt.gz", "")

        gsm_id = sample_id.split("_")[0]

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
        adata.obs["GSM"] = gsm_id

        adatas.append(adata)

    if not adatas:
        raise FileNotFoundError(
            f"No .counts.txt.gz files found in {dataset_dir}"
        )

    combined = ad.concat(
        adatas,
        join="outer",
        merge="same",
        fill_value=0,
        index_unique="-"
    )

    combined.var_names_make_unique()
    combined.obs_names_make_unique()

    return combined
