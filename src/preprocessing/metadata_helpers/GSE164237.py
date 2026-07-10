import re


def get_patient_id(sample):

    match = re.search(r"TCR-\d+", sample)

    if match:
        return match.group(0)

    return "Unknown"

def get_sample_id(sample):

    match = re.search(r"TCR-\d+-(BL|PT)", sample)

    if match:
        return match.group(0)

    return "Unknown"


def get_timepoint(sample):

    if "-BL" in sample:
        return "Baseline"

    elif "-PT" in sample:
        return "Post-Treatment"

    return "Unknown"


def get_response(sample):

    if "_R_" in sample:
        return "R"

    elif "_NR_" in sample:
        return "NR"

    return "Unknown"


def add_metadata(adata):

    adata.obs["Cancer type"] = "Melanoma"
    adata.obs["Tissue"] = "PBMC"

    adata.obs["Patient ID"] = (
        adata.obs["sample_id"]
        .apply(get_patient_id)
    )

    adata.obs["Response"] = (
        adata.obs["sample_id"]
        .apply(get_response)
    )

    adata.obs["Timepoint"] = (
        adata.obs["sample_id"]
        .apply(get_timepoint)
    )

    adata.obs["Sample ID"] = (
        adata.obs["sample_id"]
        .apply(get_sample_id)
    )

    adata.obs["Treatment"] = "Pembrolizumab"

    return adata