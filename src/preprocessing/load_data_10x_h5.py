import scanpy as sc
import anndata as ad


def load_single_10x_h5(file_path):

    # ============================================================
    # Load 10x H5 Gene Expression Matrix
    # ============================================================

    adata = sc.read_10x_h5(
        file_path
    )

    adata.var_names_make_unique()

    # Add basic metadata

    adata.obs["sample_path"] = str(file_path)
    adata.obs["sample_id"] = file_path.stem

    print(
        f"Dataset dimensions: "
        f"{adata.n_obs} cells x {adata.n_vars} genes"
    )

    return adata


def find_10x_h5_samples(input_dir):

    """
    Find all .h5 files recursively.
    """

    files = sorted(
        input_dir.rglob("*.h5")
    )

    return files


def read_multiple_10x_h5_samples(input_dir):

    """
    Load all H5 samples and combine them.
    """

    adatas = []

    for file in find_10x_h5_samples(input_dir):

        print(f"Loading: {file}")

        adata = load_single_10x_h5(file)

        adatas.append(adata)
        break

    combined = ad.concat(
        adatas,
        axis=0,
        join="outer",
        fill_value=0,
        label="sample_id_concat",
        keys=[
            a.obs["sample_id"].iloc[0]
            for a in adatas
        ],
        index_unique=None
    )

    return combined