
import os
import anndata as ad
import pandas as pd


def load_gse166181(dataset_dir):
    file = "GSE166181_Raw_UMI_CountMatrix_unnormalized_.tsv"
    path = os.path.join(dataset_dir, file)

    # Rows are genes and columns are cells
    df = pd.read_csv(
        path,
        sep="\t",
        index_col=0,
        low_memory=False
    )

    # Convert expression values to numeric
    df = df.apply(pd.to_numeric, errors="raise")

    # Transpose so cells are rows and genes are columns
    df = df.T

    adata = ad.AnnData(
        X=df.to_numpy(),
        obs=pd.DataFrame(index=df.index.astype(str)),
        var=pd.DataFrame(index=df.columns.astype(str))
    )

    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    adata.obs["sample_id"] = adata.obs_names.astype(str)
    adata.obs["source_count_file"] = path
    adata.obs["sample_path"] = path

    adata.uns["input_format"] = "tabular"
    adata.uns["input_count_files"] = [path]

    print(f"Loaded shape: {adata.shape}")
    print(f"Expression dtype: {adata.X.dtype}")

    return adata
