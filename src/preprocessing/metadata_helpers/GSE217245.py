# ============================================================
# Metadata Functions
# GSE217245 - Multiple Myeloma Anti-BCMA TCE
# ============================================================

def extract_gsm_id(sample):
    return str(sample).split("_")[0]


def extract_patient_id(sample):
    return str(sample).split("_")[1]


def extract_sample_id(sample):
    return str(sample).split("_")[-1]


def extract_response(sample):
    sample = str(sample)

    if "non_responder" in sample:
        return "NR"
    elif "responder" in sample:
        return "R"

    return "Unknown"


def extract_timepoint(sample):
    sample = str(sample)

    if "_late_relapse_" in sample:
        return "Late Relapse"
    elif "_late_" in sample:
        return "Late"
    elif "_pre_" in sample:
        return "Pre"

    return "Unknown"


def add_metadata(adata):
    """
    Add GSE217245 metadata columns to AnnData.
    """

    adata.obs["Cancer type"] = "Multiple Myeloma"
    adata.obs["Tissue"] = "Bone Marrow"

    adata.obs["Patient ID"] = (
        adata.obs["sample_id"]
        .apply(extract_patient_id)
    )

    adata.obs["Response"] = (
        adata.obs["sample_id"]
        .apply(extract_response)
    )

    adata.obs["Timepoint"] = (
        adata.obs["sample_id"]
        .apply(extract_timepoint)
    )

    adata.obs["Sample ID"] = (
        adata.obs["sample_id"]
        .apply(extract_sample_id)
    )

    adata.obs["GSM ID"] = (
        adata.obs["sample_id"]
        .apply(extract_gsm_id)
    )

    adata.obs["Treatment"] = "Anti-BCMA TCE"

    return adata