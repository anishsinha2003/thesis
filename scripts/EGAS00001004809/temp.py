import scanpy as sc
import pandas as pd
from scipy.io import mmread

counts = mmread("/Users/anishsinha/Desktop/thesis/preprocessing/datasets_to_preprocess/original_raw_datasets/EGAS00001004809/cohort2.mtx").T.tocsr()

genes = pd.read_csv(
    "/Users/anishsinha/Desktop/thesis/preprocessing/datasets_to_preprocess/original_raw_datasets/EGAS00001004809/cohort2_genes.tsv",
    header=None
)[0]

barcodes = pd.read_csv(
    "/Users/anishsinha/Desktop/thesis/preprocessing/datasets_to_preprocess/original_raw_datasets/EGAS00001004809/cohort2_barcodes.tsv",
    header=None
)[0]

meta = pd.read_csv(
    "/Users/anishsinha/Desktop/thesis/preprocessing/datasets_to_preprocess/original_raw_datasets/EGAS00001004809/1871-BIOKEY_metaData_cohort2_web.csv",
    index_col=0
)

adata = sc.AnnData(X=counts)

adata.obs_names = barcodes
adata.var_names = genes

adata.obs = meta.loc[adata.obs_names]

adata.write_h5ad("cohort2.h5ad")