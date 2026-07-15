import os

import anndata as ad
import pandas as pd
from scipy.sparse import csr_matrix


def load_single_rawcounts(file_path):

    # df = pd.read_csv(
    #     file_path,
    #     sep="\t"
    # )
    df = pd.read_csv(file_path, sep="\t", nrows=5)

    print(df.head())
    print(df.columns[:20])
    gene_symbols = df["Gene Symbol"]

    count_matrix = df.iloc[:, 3:]

    count_matrix = count_matrix.apply(
    pd.to_numeric,
    errors="coerce"
)

    count_matrix = count_matrix.fillna(0)

    adata = ad.AnnData(
        X=csr_matrix(count_matrix.T)
    )

    adata.var_names = gene_symbols.astype(str)
    adata.obs_names = count_matrix.columns.astype(str)

    adata.var["Gene_ID_1"] = df["Gene ID_1"].values
    adata.var["Gene_ID_2"] = df["Gene ID_2"].values

    adata.var_names_make_unique()

    return adata


def load_gse130157(input_dir):

    adatas = []

    for file_name in os.listdir(input_dir):

        if file_name.endswith(".RawCounts.txt"):
            print("PASS3")
            file_path = os.path.join(
                input_dir,
                file_name
            )

            print(f"Loading {file_name}")

            adata = load_single_rawcounts(
                file_path
            )
            print("PASS1")
            adatas.append(
                adata
            )
            print("PASS2")

    if len(adatas) == 0:
        raise ValueError(
            f"No RawCounts files found in {input_dir}"
        )

    adata = ad.concat(
        adatas,
        join="outer",
        merge="same"
    )

    return adata