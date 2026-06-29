# save an AnnData object to an h5ad file

import anndata as ad
import scanpy as sc

def save_h5ad_file(adata, file_path):
    # CHANGED 
    # check the obs data types and manually convert to strings
    for col in adata.obs.columns:
        if adata.obs[col].dtype == 'object':
            adata.obs[col] = adata.obs[col].astype(str)

    print(f"Saving adata to {file_path}")
    adata.write_h5ad(file_path)
    print(f"Saved adata to {file_path}")
    return