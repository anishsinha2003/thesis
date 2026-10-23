
import pandas as pd

from src.preprocessing.config import DATASETS_TO_PREPROCESS


METADATA_FILE = (
    DATASETS_TO_PREPROCESS
    / "GSE207422"
    / "GSE207422_NSCLC_scRNAseq_metadata.xlsx"
)


def get_patient_id(sample, metadata):
    row = metadata.loc[metadata["Sample"] == sample]

    if not row.empty:
        return str(row.iloc[0]["Patient"])

    return "Unknown"


def get_sample_metadata(sample, metadata):

    row = metadata.loc[metadata["Sample"] == sample]

    if row.empty:
        return {
            "Patient ID": "Unknown",
            "Sex": "Unknown",
            "Age": "Unknown",
            "Clinical Stage": "Unknown",
            "Cancer type": "NSCLC",
            "Tissue": "Unknown",
            "PD1 Antibody": "Unknown",
            "Chemotherapy": "Unknown",
            "Pathologic Response": "Unknown",
            "Residual Tumor": "Unknown",
            "RECIST": "Unknown",
            "Response": "Unknown",
        }

    row = row.iloc[0]

    recist = str(row["RECIST"]).strip().upper()

    # Binary response based on RECIST.
    if recist in {"CR", "PR"}:
        response = "R"
    elif recist in {"SD", "PD"}:
        response = "NR"
    else:
        response = "Unknown"

    return {
        "Patient ID": str(row["Patient"]),
        "Sex": row["Sex"],
        "Age": row["Age"],
        "Clinical Stage": row["Clinical Stage"],
        "Cancer type": "NSCLC",
        "Tissue": "Unknown",
        "PD1 Antibody": row["PD1 Antibody"],
        "Chemotherapy": row["Chemotherapy"],
        "Pathologic Response": row["Pathologic Response"],
        "Residual Tumor": row["Residual Tumor"],
        "RECIST": recist,
        "Response": response,
    }


def add_metadata(adata):

    metadata = pd.read_excel(METADATA_FILE)

    # Clean the sample identifiers in the metadata.
    metadata["Sample"] = (
        metadata["Sample"].astype(str).str.strip()
    )

    # Ensure each cell has a sample_id matching the Excel
    # metadata's Sample column.
    if "sample_id" not in adata.obs.columns:
        raise KeyError(
            "adata.obs must contain 'sample_id' with the "
            "original sample names, such as BD_immune01."
        )

    sample_ids = adata.obs["sample_id"].astype(str).str.strip()

    metadata_lookup = metadata.set_index("Sample")

    for column in [
        "Patient",
        "Sex",
        "Age",
        "Clinical Stage",
        "Pathology",
        "PD1 Antibody",
        "Chemotherapy",
        "Pathologic Response",
        "Residual Tumor",
        "RECIST",
    ]:
        if column in metadata_lookup.columns:
            adata.obs[column] = (
                sample_ids.map(metadata_lookup[column])
                .fillna("Unknown")
            )

    adata.obs["Patient ID"] = (
        adata.obs["Patient"].fillna("Unknown").astype(str)
    )

    adata.obs["Sample ID"] = sample_ids

    adata.obs["Cancer type"] = "NSCLC"
    adata.obs["Tissue"] = "Unknown"

    adata.obs["Response"] = (
        adata.obs["RECIST"]
        .astype(str)
        .str.strip()
        .str.upper()
        .map({
            "CR": "R",
            "PR": "R",
            "SD": "NR",
            "PD": "NR",
        })
        .fillna("Unknown")
    )

    return adata
