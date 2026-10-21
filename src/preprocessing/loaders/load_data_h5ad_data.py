
# Read h5ad files into AnnData objects

import anndata as ad
import scanpy as sc


def read_single_h5ad_file(h5ad_file):
    adata = sc.read_h5ad(h5ad_file)

    # Keep track of the file and sample
    adata.obs["sample_path"] = str(h5ad_file)
    adata.obs["sample_id"] = h5ad_file.stem

    return adata


def read_multiple_h5ad_files(input_dir):
    # Find all h5ad files in the input directory
    h5ad_files = sorted(input_dir.glob("*.h5ad"))

    if not h5ad_files:
        raise FileNotFoundError(
            f"No h5ad files found in {input_dir}"
        )

    # Read each file
    adatas = [
        read_single_h5ad_file(file)
        for file in h5ad_files
    ]

    # Combine all samples
    combined = ad.concat(
        adatas,
        axis=0,
        join="outer",
        fill_value=0,
        label="sample_id_concat",
        keys=[adata.obs["sample_id"].iloc[0] for adata in adatas],
        index_unique=None,
    )

    return combined
