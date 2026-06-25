# src/preprocessing/preprocess_dataset.py

from pathlib import Path

from src.preprocessing.load_data_h5ad_data import read_multiple_h5ad_files
from src.preprocessing.load_data_tabular import load_tabular_folder
from src.preprocessing.load_data_10x_h5 import read_multiple_10x_h5_samples
from src.preprocessing.load_data_10x_mtx import read_multiple_10x_samples
from src.preprocessing.qc import run_qc
from src.preprocessing.normalise import normalise
from src.preprocessing.run_celltypist import run_celltypist
from src.preprocessing.align_gene_to_hg38 import align_gene_to_hg38, preprocess_gtf_file_hg38
from src.preprocessing.save_h5ad_file import save_h5ad_file
# from src.preprocessing.preprocess_dataset import preprocess_dataset

# TODO:
# create this function later
# from src.preprocessing.filter_immune_cells import filter_immune_cells

def validate_processed_adata(adata, hg38_gene_df):

    print("=" * 60)
    print("ADATA VALIDATION")
    print("=" * 60)

    # ------------------------------------------------
    # Basic dimensions
    # ------------------------------------------------

    print(f"Cells : {adata.n_obs}")
    print(f"Genes : {adata.n_vars}")

    # ------------------------------------------------
    # Duplicate cells
    # ------------------------------------------------

    print(
        f"Duplicate cell names: "
        f"{adata.obs_names.duplicated().sum()}"
    )

    # ------------------------------------------------
    # Duplicate genes
    # ------------------------------------------------

    print(
        f"Duplicate gene names: "
        f"{adata.var_names.duplicated().sum()}"
    )

    # ------------------------------------------------
    # Layers
    # ------------------------------------------------

    print("\nLayers:")
    print(list(adata.layers.keys()))

    if "counts" in adata.layers:
        print("✓ counts layer exists")
    else:
        print("✗ counts layer missing")

    # ------------------------------------------------
    # QC columns
    # ------------------------------------------------

    qc_cols = [
        "n_genes_by_counts",
        "total_counts",
    ]

    print("\nQC Columns")

    for col in qc_cols:
        if col in adata.obs.columns:
            print(f"✓ {col}")
        else:
            print(f"✗ {col}")

    # ------------------------------------------------
    # CellTypist
    # ------------------------------------------------

    print("\nCell Type Annotation")

    if "predicted_labels" in adata.obs.columns:
        print("✓ CellTypist labels found")
    else:
        print("✗ CellTypist labels missing")

    # ------------------------------------------------
    # Normalisation
    # ------------------------------------------------

    print("\nNormalisation")

    if "log1p" in adata.uns:
        print("✓ log1p found")
    else:
        print("✗ log1p not found")

    # ------------------------------------------------
    # HG38 overlap
    # ------------------------------------------------

    hg38_overlap = adata.var_names.isin(
        hg38_gene_df["gene_symbol"]
    )

    pct = (
        hg38_overlap.sum()
        / adata.n_vars
        * 100
    )

    print("\nHG38 Alignment")

    print(
        f"{hg38_overlap.sum()} / "
        f"{adata.n_vars}"
    )

    print(
        f"Overlap = {pct:.2f}%"
    )

    # ------------------------------------------------
    # Missing values
    # ------------------------------------------------

    print("\nMissing Values")

    print(
        f"obs NaNs: "
        f"{adata.obs.isna().sum().sum()}"
    )

    print(
        f"var NaNs: "
        f"{adata.var.isna().sum().sum()}"
    )

    print("=" * 60)


def preprocess_dataset(
    input_dir,
    output_file,
    data_format,
    prefix=""
):
    """
    Complete preprocessing pipeline for one dataset.
    """

    print("=" * 60)
    print("LOADING DATA")
    print("=" * 60)

    # ============================================================
    # Load Dataset
    # ============================================================

    if data_format == "10x":
        adata = read_multiple_10x_samples(
            Path(input_dir),
            prefix=prefix
        )

    elif data_format == "h5ad":
        adata = read_multiple_h5ad_files(
            Path(input_dir)
        )

    elif data_format == "tabular":
        adata = load_tabular_folder(
            Path(input_dir)
        )
    elif data_format == "h5":
        adata = read_multiple_10x_h5_samples(
            Path(input_dir)
        )

    else:
        raise ValueError(
            f"Unsupported format: {data_format}"
        )

    print(adata)

    # ============================================================
    # Normalisation
    # ============================================================

    print("\nRunning Normalisation...")
    adata = normalise(adata)

    print(adata)

    # ============================================================
    # QC
    # ============================================================

    print("\nRunning QC...")
    adata = run_qc(adata)

    print(adata)

    # ============================================================
    # Cell Type Annotation
    # ============================================================

    print("\nRunning CellTypist...")
    adata = run_celltypist(adata)

    print(adata)

    # ============================================================
    # Keep Immune Cells
    # ============================================================

    print("\nFiltering Immune Cells...")

    # TODO:
    # adata = filter_immune_cells(adata)

    print(adata)

    # ============================================================
    # Gene Alignment
    # ============================================================

    print("\nAligning Genes To HG38...")
    adata = align_gene_to_hg38(adata)

    # ============================================================
    # Validate Processed Data
    # ============================================================
    print("\nValidating Processed AnnData...")
    hg38_gene_df = preprocess_gtf_file_hg38()

    validate_processed_adata(
        adata,
        hg38_gene_df
    )

    # ============================================================
    # Save Dataset
    # ============================================================

    save_h5ad_file(adata, output_file)

    return adata


preprocess_dataset(
    input_dir="/Users/anishsinha/Desktop/thesis/preprocessing/datasets_to_preprocess/GSE217245",
    output_file="/Users/anishsinha/Desktop/thesis/preprocessing/preprocessed_data/GSE217245",
    data_format="10x",
    prefix=""
)
