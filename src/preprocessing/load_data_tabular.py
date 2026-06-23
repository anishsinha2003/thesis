# read tabular data file into an AnnData object (tsv, csv, txt)
# this can link to the R code that was used to generate the h5ad files

import anndata as ad
import scanpy as sc
import pandas as pd

def guess_sep(path) -> str:
    path = str(path).lower()
    if path.endswith(".csv") or path.endswith(".csv.gz"):
        return ","
    return "\t"

def read_tabular_file(file_path):
    sep = guess_sep(file_path)
    df = pd.read_csv(file_path, sep=sep, index_col=0)
    adata = ad.AnnData(X=df.values)
    adata.obs_names = df.index.astype(str)
    adata.var_names = df.columns.astype(str)
    adata.var_names_make_unique()
    adata.obs["source_count_file"] = str(file_path)
    adata.obs["sample_id"] = file_path.stem.replace(".tsv", "").replace(".csv", "").replace(".txt", "")
    adata.obs["sample_path"] = file_path
    adata.uns["input_format"] = "tabular"
    adata.uns["input_count_files"] = [str(file_path)]
    return adata

def find_count_tables(input_dir):
    """
    Find tabular count files in a folder.
    """

    patterns = ["*.csv", "*.tsv", "*.txt", "*.csv.gz", "*.tsv.gz", "*.txt.gz"]
    files = []

    for pattern in patterns:
        files.extend(input_dir.rglob(pattern))

    files = sorted(files)

    if len(files) == 0:
        raise FileNotFoundError(f"No tabular count files found in: {input_dir}")

    return files


def load_tabular_folder(input_dir):
    # read all tabular files in the input directory
    tabular_files = find_count_tables(input_dir)
    adatas = [read_tabular_file(f) for f in tabular_files]
    return ad.concat(adatas, axis=0, join="outer", fill_value=0, label="sample_id_concat", keys=[a.obs["sample_id"].iloc[0] for a in adatas], index_unique=None)