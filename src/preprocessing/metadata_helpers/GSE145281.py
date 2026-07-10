import re


def get_patient_id(sample):

    match = re.search(r"(R\d+|NR\d+)", sample)

    if match:
        return match.group(1)

    return "Unknown"


def get_sample_id(sample):

    match = re.search(r"GSM\d+", sample)

    if match:
        return match.group(0)

    return "Unknown"


def get_response(sample):

    if "_NR" in sample:
        return "NR"

    elif "_R" in sample:
        return "R"

    return "Unknown"


def add_metadata(adata):

    adata.obs["Cancer type"] = (
        "Metastatic Urothelial Carcinoma"
    )

    adata.obs["Tissue"] = "PBMC"

    adata.obs["Patient ID"] = (
        adata.obs["sample_id"]
        .apply(get_patient_id)
    )

    adata.obs["Response"] = (
        adata.obs["sample_id"]
        .apply(get_response)
    )

    adata.obs["Sample ID"] = (
        adata.obs["sample_id"]
        .apply(get_sample_id)
    )

    adata.obs["Treatment"] = "Atezolizumab"

    return adata