# save an AnnData object to an h5ad file

import anndata as ad
import scanpy as sc

def save_h5ad_file(adata, file_path):
    print(f"Saving adata to {file_path}")
    adata.write_h5ad(file_path)
    print(f"Saved adata to {file_path}")
    return