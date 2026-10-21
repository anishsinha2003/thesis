
import re


def get_patient_id(sample):

    match = re.search(r"Patient(\d+)", sample)

    if match:
        return f"Patient{match.group(1)}"

    return "Unknown"


def get_sample_id(sample):

    match = re.search(r"GSM\d+", sample)

    if match:
        return match.group(0)

    return "Unknown"


def get_response(sample):

    # The primary tumour progressed during PD-1 treatment.
    # Do not assume this represents the response of the whole patient.

    if "Patient1" in sample:
        return "Unknown"

    return "Unknown"


def add_metadata(adata):

    adata.obs["Cancer type"] = (
        "MSI-H Colorectal Cancer"
    )

    adata.obs["Tissue"] = "Tumor"

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

    adata.obs["Treatment"] = "PD-1 Inhibitor"

    adata.obs["Time Point"] = "Unknown"

    return adata
