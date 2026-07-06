# ============================================================
# Import Required Libraries
# ============================================================

# scanpy:
# Main framework used for single-cell RNA-seq analysis.
# Provides tools for preprocessing, quality control,
# dimensionality reduction, graph construction, and visualisation.

import scanpy as sc

# pandas:
# Used for reading and handling tabular files such as
# gene lists and barcode metadata.

import pandas as pd

import anndata as ad

# mmread:
# Used to load Matrix Market (.mtx) sparse matrix files,
# which store the raw gene expression counts.

import scanpy as sc

def load_single_10x_sample(sample_path, prefix):

    # ============================================================
    # Load 10x Mtx Gene Expression Matrix
    # ============================================================

    # The matrix file contains raw gene expression counts.
    # Rows originally represent genes and columns represent cells.
    # The matrix is transposed (.T) so that:
    #
    # rows   = cells
    # columns = genes
    #
    # which is the format expected by Scanpy.

    adata = sc.read_10x_mtx(
        sample_path,
        var_names="gene_symbols",
        make_unique=True,
        prefix=prefix,
    )

    # add basic metadata column to the AnnData object
    adata.obs["sample_path"] = str(sample_path)
    adata.obs["sample_id"] = sample_path.name

    # Print dataset dimensions:
    print(f"Dataset dimensions: {adata.n_obs} cells x {adata.n_vars} genes")
    # n_obs  = number of cells
    # n_vars = number of genes
    return adata

def is_10x_mtx_folder(folder, prefix):
    """
    Check whether a folder contains standard 10x mtx files.
    Supports optional prefix used by scanpy.read_10x_mtx.
    """

    files = {p.name for p in folder.iterdir() if p.is_file()}

    matrix_names = {f"{prefix}matrix.mtx", f"{prefix}matrix.mtx.gz"}
    barcode_names = {f"{prefix}barcodes.tsv", f"{prefix}barcodes.tsv.gz"}
    feature_names = {
        f"{prefix}features.tsv", f"{prefix}features.tsv.gz",
        f"{prefix}genes.tsv", f"{prefix}genes.tsv.gz",
    }

    return (
        len(files & matrix_names) > 0
        and len(files & barcode_names) > 0
        and len(files & feature_names) > 0
    )

def find_10x_mtx_samples(input_dir, prefix):
    """
    Recursively find sample folders containing standard 10x mtx files.
    """

    folders = []
    if is_10x_mtx_folder(input_dir, prefix=prefix):
        folders.append(input_dir)

    for folder in input_dir.rglob("*"):
        if folder.is_dir() and is_10x_mtx_folder(folder, prefix=prefix):
            folders.append(folder)

    folders = sorted(set(folders))
    return folders

def read_multiple_10x_samples(input_dir, prefix):
    # parse the sample dir, each dir has multiple 10x mtx directories
    # read each 10x mtx directory into an AnnData object
    # concatenate the AnnData objects into a single AnnData object
    # return the AnnData object

    adatas = []

    for folder in find_10x_mtx_samples(input_dir, prefix):
        adata = load_single_10x_sample(folder, prefix)
        adatas.append(adata)
        # add the path to the sample to the AnnData object

    return ad.concat(adatas, axis=0, join="outer", fill_value=0, label="sample_id_concat", keys=[a.obs["sample_path"].iloc[0] for a in adatas], index_unique=None)