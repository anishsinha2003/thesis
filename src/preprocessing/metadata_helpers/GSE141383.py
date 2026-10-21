
# ============================================================
# GSE141383 Metadata
# Glioma / Glioblastoma
# ============================================================

import pandas as pd


# ============================================================
# GEO Sample Metadata
# ============================================================

SAMPLE_METADATA = {

    "GSM4202144": {
        "sample": "PJ052",
    },

    "GSM4202145": {
        "sample": "PJ053",
    },

    "GSM4202146": {
        "sample": "PW016-703",
    },

    "GSM4202147": {
        "sample": "PW017-703",
    },

    "GSM4202148": {
        "sample": "PW032-706",
    },

    "GSM4202149": {
        "sample": "PW032-710",
    },

    "GSM5975476": {
        "sample": "PW035-710",
    },

    "GSM5975477": {
        "sample": "PDC001",
    },

    "GSM5975478": {
        "sample": "PW039-705",
    },
}


# ============================================================
# Helper Functions
# ============================================================

def get_gsm(sample):

    import re

    match = re.search(r"GSM\d+", str(sample))

    if match:
        return match.group(0)

    return "Unknown"


def get_patient_id(sample):

    gsm = get_gsm(sample)

    return SAMPLE_METADATA.get(
        gsm, {}
    ).get("sample", "Unknown")


def get_sample_id(sample):

    return str(sample)


# ============================================================
# Metadata Alignment
# ============================================================

def add_metadata(adata):

    adata.obs["Cancer type"] = (
        "Glioma"
    )

    adata.obs["Tissue"] = (
        "Glioma surgical specimen"
    )

    adata.obs["Patient ID"] = (
        adata.obs["sample_id"]
        .apply(get_patient_id)
    )

    adata.obs["Sample ID"] = (
        adata.obs["sample_id"]
        .apply(get_sample_id)
    )

    # These labels are not assigned here because they
    # should not be inferred from the sample filenames.
    adata.obs["Response"] = "Unknown"
    adata.obs["Treatment"] = "Unknown"

    return adata
