import pandas as pd
import re

from ..config import DATASETS_TO_PREPROCESS

# ============================================================
# GSE212217 Metadata
# Endometrial Cancer
# Pembrolizumab
# ============================================================

# Load metadata file
metadata = pd.read_csv(
    DATASETS_TO_PREPROCESS / "GSE212217/GSE212217_seurat_scRNAseq_metadata.txt",
    sep="\t"
)



# Create patient -> clinical mapping
response_lookup = (
    metadata[["patient", "clinical"]]
    .drop_duplicates()
    .set_index("patient")["clinical"]
    .to_dict()
)

def extract_sample_id(sample):
    sample = str(sample)

    match = re.search(r"(PEM\d+C\d+)", sample)

    if match:
        return match.group(1)

    return sample

def extract_patient_id(sample):

    sample = str(sample)

    match = re.search(r"PEM(\d+)", sample)

    if match:
        return int(match.group(1))

    return None


def get_response(sample):

    patient = extract_patient_id(sample)

    if patient is None:
        return "Unknown"

    response = response_lookup.get(patient, "Unknown")

    if response in ["epiR", "mutR"]:
        return "R"
    elif response == "NR":
        return "NR"

    return "Unknown"


def get_timepoint(sample):

    sample = str(sample)

    if "C1" in sample:
        return "Pre"

    elif "C3" in sample:
        return "On Treatment"

    elif "C5" in sample:
        return "On Treatment"

    return "Unknown"


def add_metadata(adata):

    adata.obs["Cancer type"] = "Endometrial Carcinoma"

    adata.obs["Tissue"] = "PBMC"

    adata.obs["Treatment"] = "Pembrolizumab"

    adata.obs["Treatment Class"] = "Anti-PD-1"

    adata.obs["Sample ID"] = (
        adata.obs["sample_id"]
        .apply(extract_sample_id)
    )


    adata.obs["Patient ID"] = (
        adata.obs["sample_id"]
        .apply(extract_patient_id)
    )

    adata.obs["Response"] = (
        adata.obs["sample_id"]
        .apply(get_response)
    )

    adata.obs["Timepoint"] = (
        adata.obs["sample_id"]
        .apply(get_timepoint)
    )

    return adata