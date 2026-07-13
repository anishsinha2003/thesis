import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix


def load_gse130157(file_path):

    df = pd.read_csv(
        file_path,
        sep="\t"
    )

    gene_symbols = df["Gene_Symbol"]

    count_matrix = df.iloc[:, 3:]

    adata = ad.AnnData(
        X=csr_matrix(count_matrix.T)
    )

    adata.var_names = gene_symbols.astype(str)
    adata.obs_names = count_matrix.columns.astype(str)

    adata.var["Gene_ID_1"] = (
        df["Gene_ID_1"].values
    )

    adata.var["Gene_ID_2"] = (
        df["Gene_ID_2"].values
    )

    adata.var_names_make_unique()

    return adata