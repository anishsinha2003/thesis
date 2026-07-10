import re

PATIENT_METADATA = {
    "P01": {
        "response": "NR",
        "therapy": "NAC",
        "pathology": "SCC",
        "age": 68,
        "sex": "Male",
    },
    "P02": {
        "response": "NR",
        "therapy": "NAPC",
        "pathology": "SCC",
        "age": 69,
        "sex": "Male",
    },
    "P03": {
        "response": "R",
        "therapy": "NAPC",
        "pathology": "SCC",
        "age": 51,
        "sex": "Female",
    },
    "P04": {
        "response": "NR",
        "therapy": "NAPC",
        "pathology": "AD",
        "age": 66,
        "sex": "Female",
    },
    "P05": {
        "response": "R",
        "therapy": "NAPC",
        "pathology": "SCC",
        "age": 67,
        "sex": "Male",
    },
    "P06": {
        "response": "R",
        "therapy": "NAPC",
        "pathology": "SCC",
        "age": 62,
        "sex": "Male",
    },
    "P07": {
        "response": "NR",
        "therapy": "NAPC",
        "pathology": "AD",
        "age": 46,
        "sex": "Male",
    },
}


import re

def get_patient_id(sample):

    match = re.search(r"P\d+", sample)

    if match:
        return match.group(0)

    return "Unknown"


def get_response(sample):

    patient = get_patient_id(sample)

    return (
        PATIENT_METADATA
        .get(patient, {})
        .get("response", "Unknown")
    )


def get_therapy(sample):

    patient = get_patient_id(sample)

    return (
        PATIENT_METADATA
        .get(patient, {})
        .get("therapy", "Unknown")
    )


def get_pathology(sample):

    patient = get_patient_id(sample)

    return (
        PATIENT_METADATA
        .get(patient, {})
        .get("pathology", "Unknown")
    )


def get_age(sample):

    patient = get_patient_id(sample)

    return (
        PATIENT_METADATA
        .get(patient, {})
        .get("age", None)
    )


def get_sex(sample):

    patient = get_patient_id(sample)

    return (
        PATIENT_METADATA
        .get(patient, {})
        .get("sex", "Unknown")
    )


def add_metadata(adata):

    adata.obs["Cancer type"] = "Non-Small Cell Lung Cancer"
    adata.obs["Tissue"] = "Tumor"

    adata.obs["Patient ID"] = (
        adata.obs["sample_id"]
        .apply(get_patient_id)
    )

    adata.obs["Sample ID"] = (
        adata.obs["sample_id"]
    )

    adata.obs["Response"] = (
        adata.obs["sample_id"]
        .apply(get_response)
    )

    adata.obs["Treatment"] = (
        adata.obs["sample_id"]
        .apply(get_therapy)
    )

    adata.obs["Pathology"] = (
        adata.obs["sample_id"]
        .apply(get_pathology)
    )

    adata.obs["Age"] = (
        adata.obs["sample_id"]
        .apply(get_age)
    )

    adata.obs["Sex"] = (
        adata.obs["sample_id"]
        .apply(get_sex)
    )

    return adata