# read h5ad file into an AnnData object
# this can link to the R code that was used to generate the h5ad files

import anndata as ad
import scanpy as sc

def read_single_h5ad_file(h5ad_file):
    adata = sc.read_h5ad(h5ad_file)
    adata.obs["sample_path"] = h5ad_file
    adata.obs["sample_id"] = h5ad_file.name
    return adata
    
def read_multiple_h5ad_files(input_dir):
    # read all h5ad files in the input directory
    h5ad_files = [f for f in input_dir.glob("*.h5ad")]
    adatas = [read_single_h5ad_file(f) for f in h5ad_files]
    return ad.concat(adatas, axis=0, join="outer", fill_value=0, label="sample_id_concat", keys=[a.obs["sample_id"].iloc[0] for a in adatas], index_unique=None)