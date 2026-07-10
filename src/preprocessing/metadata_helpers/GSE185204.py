import re


def get_patient_id(sample):

    match = re.search(
        r"(MSK\d+)",
        sample
    )

    if match:
        return match.group(1)

    return "Unknown"


def get_sample_id(sample):

    sample = sample.replace(".h5", "")

    return sample


def get_tissue(sample):

    if "_LN_" in sample:
        return "Lymph Node"

    if "_Normal_" in sample:
        return "Adjacent Normal"

    return "Tumor"


def get_response(sample):

    patient = get_patient_id(sample)

    response_map = {
        "MSK1263": "NR",
        "MSK1302": "R",
        "MSK1344": "Intermediate",
    }

    return response_map.get(
        patient,
        "Unknown"
    )


def get_treatment(sample):

    return "Anti-PD-1"


def add_metadata(adata):

    adata.obs["Cancer type"] = (
        "Non-Small Cell Lung Cancer"
    )

    adata.obs["Patient ID"] = (
        adata.obs["sample_id"]
        .apply(get_patient_id)
    )

    adata.obs["Sample ID"] = (
        adata.obs["sample_id"]
        .apply(get_sample_id)
    )

    adata.obs["Response"] = (
        adata.obs["sample_id"]
        .apply(get_response)
    )

    adata.obs["Tissue"] = (
        adata.obs["sample_id"]
        .apply(get_tissue)
    )

    adata.obs["Treatment"] = (
        adata.obs["sample_id"]
        .apply(get_treatment)
    )

    adata.obs["Treatment Class"] = (
        "Anti-PD-1"
    )

    return adata