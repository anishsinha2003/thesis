# src/preprocessing/preprocess_dataset.py

from pathlib import Path

from src.preprocessing.dataset_config import DATASETS
from src.preprocessing.config import GTF_FILE, GFF3_PATH
from src.preprocessing.loaders.load_data_h5ad_data import read_multiple_h5ad_files
from src.preprocessing.loaders.load_data_tabular import load_tabular_folder
from src.preprocessing.loaders.load_data_10x_h5 import read_multiple_10x_h5_samples
from src.preprocessing.loaders.load_data_10x_mtx import read_multiple_10x_samples
from src.preprocessing.loaders.dataset_specfic_loaders.GSE169246 import load_gse169246
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
    # Samples
    # ------------------------------------------------

    print("\nSamples")

    if "Sample ID" in adata.obs.columns:

        sample_counts = adata.obs["Sample ID"].value_counts()

        print(f"Number of unique samples: {len(sample_counts)}")

        if len(sample_counts) > 20:

            for sample, count in sample_counts.iloc[:20].items():
                print(f"{sample}: {count} cells")

            print("...")

        else:

            for sample, count in sample_counts.items():
                print(f"{sample}: {count} cells")

    else:
        print("✗ Sample ID column missing")

    required_cols = ["Patient ID", "Response"]

    if all(col in adata.obs.columns for col in required_cols):

        patient_summary = (
            adata.obs[
                ["Patient ID", "Response"]
            ]
            .drop_duplicates()
            .copy()
        )

        patient_summary["Cell Count"] = (
            patient_summary["Patient ID"]
            .map(
                adata.obs.groupby("Patient ID")
                .size()
            )
        )

        patient_summary = patient_summary.sort_values(
            "Patient ID"
        )

        print(f"Number of unique patients: {len(patient_summary)}")
        print()

        print(
            patient_summary[
                [
                    "Patient ID",
                    "Response",
                    "Cell Count",
                ]
            ].to_string(index=False)
        )

    else:

        missing = [
            col for col in required_cols
            if col not in adata.obs.columns
        ]

        print(
            f"✗ Missing columns: {missing}"
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

def validate_required_metadata(adata):
    """
    Check that required metadata columns exist in adata.obs.
    Raises an error if any are missing.
    """

    required_columns = [
        "Cancer type",
        "Tissue",
        "Response",
        "Patient ID",
        "Sample ID",
        "Treatment",
    ]

    missing_columns = [
        col for col in required_columns
        if col not in adata.obs.columns
    ]

    if missing_columns:
        raise ValueError(
            "Missing required metadata columns: "
            + ", ".join(missing_columns)
        )

    print("✓ Metadata validation passed")
    print("✓ Required columns present:")
    print(required_columns)

def preprocess_dataset(
    input_dir,
    output_file,
    data_format,
    metadata_fn,
    prefix="",

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
    elif callable(data_format):
        adata = data_format(Path(input_dir))
    else:
        raise ValueError(
            f"Unsupported format: {data_format}"
        )

    print("PASSSSSSSSSSS ----------------------------------------------- ")
    print(adata)
    print(adata.obs_names[:20])
    print(sum("Post" in x for x in adata.obs_names))
    print(sum("Pre" in x for x in adata.obs_names))
    samples = sorted(set(x.split(".")[-1] for x in adata.obs_names))
    print(samples)

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
    gff3_path = GFF3_PATH
    gtf_file = GTF_FILE
    adata = align_gene_to_hg38(adata, gff3_path, gtf_file)

    # ============================================================
    # Metadata Alignment
    # ============================================================

    print("\MetaData Alignment...")
    adata = metadata_fn(adata)
    print(adata.obs.head().T)

    print("----------------------------------------- CHeCK AFTER thIS -----------------------------------------")
    print(sorted(adata.obs["Sample ID"].unique()))
    print(adata.obs["Timepoint"].value_counts())

    post_samples = adata.obs[
    adata.obs["Sample ID"].str.startswith("Post")
    ]["Sample ID"].unique()

    print(post_samples)
    print(f"Number of post samples: {len(post_samples)}")



    # ============================================================
    # Validate Processed Data
    # ============================================================
    print("\nValidating AnnData...")
    hg38_gene_df = preprocess_gtf_file_hg38(gff3_path)

    validate_processed_adata(
        adata,
        hg38_gene_df
    )

    validate_required_metadata(
        adata,
    )

    # ============================================================
    # Save Dataset
    # ============================================================

    save_h5ad_file(adata, output_file)

    return adata

dataset = DATASETS["GSE169246"]

preprocess_dataset(
    input_dir=dataset["input_dir"],
    output_file=dataset["output_file"],
    data_format=dataset["data_format"],
    prefix=dataset["prefix"],
    metadata_fn=dataset["metadata_fn"],
)