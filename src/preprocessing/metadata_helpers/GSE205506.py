# ============================================================
# GSE205506 Metadata
# dMMR/MSI-H Colorectal Cancer
# Anti-PD-1 ± Celecoxib
# ============================================================

import re


# ============================================================
# Patient-Level Clinical Metadata
# ============================================================

PATIENT_METADATA = {

    # Anti-PD-1

    "P12": {
        "response": "NR",
        "treatment": "Anti-PD-1",
        "sex": "M",
        "age": 67,
    },

    "P15": {
        "response": "R",
        "treatment": "Anti-PD-1",
        "sex": "F",
        "age": 58,
    },

    "P17": {
        "response": "R",
        "treatment": "Anti-PD-1",
        "sex": "M",
        "age": 62,
    },

    "P18": {
        "response": "NR",
        "treatment": "Anti-PD-1",
        "sex": "F",
        "age": 31,
    },

    "P23": {
        "response": "R",
        "treatment": "Anti-PD-1",
        "sex": "M",
        "age": 38,
    },

    "P27": {
        "response": "R",
        "treatment": "Anti-PD-1",
        "sex": "F",
        "age": 67,
    },

    "P28": {
        "response": "R",
        "treatment": "Anti-PD-1",
        "sex": "M",
        "age": 57,
    },

    "P29": {
        "response": "R",
        "treatment": "Anti-PD-1",
        "sex": "M",
        "age": 36,
    },

    "P30": {
        "response": "R",
        "treatment": "Anti-PD-1",
        "sex": "M",
        "age": 45,
    },

    "P31": {
        "response": "NR",
        "treatment": "Anti-PD-1",
        "sex": "M",
        "age": 53,
    },

    # Anti-PD-1 + Celecoxib

    "P11": {
        "response": "R",
        "treatment": "Anti-PD-1 + Celecoxib",
        "sex": "M",
        "age": 65,
    },

    "P14": {
        "response": "R",
        "treatment": "Anti-PD-1 + Celecoxib",
        "sex": "M",
        "age": 50,
    },

    "P19": {
        "response": "R",
        "treatment": "Anti-PD-1 + Celecoxib",
        "sex": "M",
        "age": 26,
    },

    "P21": {
        "response": "R",
        "treatment": "Anti-PD-1 + Celecoxib",
        "sex": "F",
        "age": 52,
    },

    "P24": {
        "response": "R",
        "treatment": "Anti-PD-1 + Celecoxib",
        "sex": "F",
        "age": 45,
    },

    "P25": {
        "response": "R",
        "treatment": "Anti-PD-1 + Celecoxib",
        "sex": "F",
        "age": 69,
    },

    "P26": {
        "response": "NR",
        "treatment": "Anti-PD-1 + Celecoxib",
        "sex": "F",
        "age": 45,
    },

    "P32": {
        "response": "R",
        "treatment": "Anti-PD-1 + Celecoxib",
        "sex": "M",
        "age": 33,
    },

    "P33": {
        "response": "R",
        "treatment": "Anti-PD-1 + Celecoxib",
        "sex": "M",
        "age": 37,
    },
}


# ============================================================
# GSM → Patient Mapping
# ============================================================

GSM_TO_PATIENT = {

    "GSM6213956": "P11",
    "GSM6213957": "P11",

    "GSM6213958": "P12",
    "GSM6213959": "P12",

    "GSM6213960": "P14",

    "GSM6213961": "P15",
    "GSM6213962": "P15",

    "GSM6213963": "P17",
    "GSM6213964": "P17",

    "GSM6213965": "P18",
    "GSM6213966": "P18",

    "GSM6213967": "P19",
    "GSM6213968": "P19",

    "GSM6213969": "P21",
    "GSM6213970": "P21",

    "GSM6213971": "P23",

    "GSM6213972": "P24",

    "GSM6213973": "P25",
    "GSM6213974": "P25",

    "GSM6213975": "P26",

    "GSM6213976": "P27",
    "GSM6213977": "P27",

    "GSM6213978": "P28",
    "GSM6213979": "P28",

    "GSM6213980": "P29",
    "GSM6213981": "P29",

    "GSM6213982": "P30",
    "GSM6213983": "P30",

    "GSM6213984": "P31",
    "GSM6213985": "P31",

    "GSM6213986": "P32",
    "GSM6213987": "P32",

    "GSM6213988": "P33",
    "GSM6213989": "P33",
}


# ============================================================
# Helper Functions
# ============================================================

def get_gsm(sample):

    match = re.search(
        r"(GSM\d+)",
        str(sample)
    )

    if match:
        return match.group(1)

    return "Unknown"


def get_patient_id(sample):

    gsm = get_gsm(sample)

    return GSM_TO_PATIENT.get(
        gsm,
        "Unknown"
    )


def get_response(sample):

    patient = get_patient_id(sample)

    return PATIENT_METADATA.get(
        patient,
        {}
    ).get(
        "response",
        "Unknown"
    )


def get_treatment(sample):

    patient = get_patient_id(sample)

    return PATIENT_METADATA.get(
        patient,
        {}
    ).get(
        "treatment",
        "Unknown"
    )


def get_age(sample):

    patient = get_patient_id(sample)

    return PATIENT_METADATA.get(
        patient,
        {}
    ).get(
        "age",
        None
    )


def get_sex(sample):

    patient = get_patient_id(sample)

    return PATIENT_METADATA.get(
        patient,
        {}
    ).get(
        "sex",
        "Unknown"
    )


def get_sample_id(sample):

    return str(sample)


# ============================================================
# Metadata Alignment
# ============================================================

def add_metadata(adata):

    adata.obs["Cancer type"] = (
        "Mismatch Repair Deficient Colorectal Cancer"
    )

    adata.obs["Tissue"] = (
        "Colon"
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

    adata.obs["Treatment"] = (
        adata.obs["sample_id"]
        .apply(get_treatment)
    )

    adata.obs["Treatment Class"] = (
        "Anti-PD-1"
    )

    adata.obs["Age"] = (
        adata.obs["sample_id"]
        .apply(get_age)
    )

    adata.obs["Sex"] = (
        adata.obs["sample_id"]
        .apply(get_sex)
    )

    for patient in sorted(
        adata.obs["Patient ID"].unique()
    ):

        response = PATIENT_METADATA.get(
            patient,
            {}
        ).get(
            "response",
            "Unknown"
        )

        print(
            f"{patient:<6} {response}"
        )

    return adata