# src/preprocessing/loaders/dataset_specific_loaders/GSE270235.py

from pathlib import Path
import scanpy as sc


def load_gse270235(input_dir):

    h5ad_files = sorted(
        Path(input_dir).glob("*.h5ad")
    )

    if len(h5ad_files) == 0:
        raise FileNotFoundError(
            f"No .h5ad files found in {input_dir}"
        )

    if len(h5ad_files) > 1:
        print(
            f"Found {len(h5ad_files)} h5ad files. "
            f"Using {h5ad_files[0].name}"
        )

    adata = sc.read_h5ad(
        h5ad_files[0]
    )

    adata.obs["sample_path"] = str(
        h5ad_files[0]
    )

    adata.obs["sample_id"] = (
        h5ad_files[0].stem
    )

    print(
        f"Dataset dimensions: "
        f"{adata.n_obs} cells x {adata.n_vars} genes"
    )

    return adata