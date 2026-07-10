import re


RESPONDERS = {
    "P019",
    "P010",
    "P012",
    "P007",
    "P022",
    "P011",
    "P020",
    "P008",
    "P013",
}

NON_RESPONDERS = {
    "P017",
    "P001",
    "P002",
    "P014",
    "P004",
    "P005",
    "P016",
    "P025",
    "P018",
    "P023",
    "P024",
    "P003",
}


def get_patient_id(sample):

    match = re.search(r"(P\d+)", sample)

    if match:
        return match.group(1)

    return "Unknown"


def get_sample_id(sample):

    match = re.search(r"(Pre_P\d+_[bt])", sample)

    if match:
        return match.group(1)

    return "Unknown"


def get_timepoint(sample):

    if "Pre_" in sample:
        return "Baseline"

    return "Unknown"


def get_response(sample):

    patient = get_patient_id(sample)

    if patient in RESPONDERS:
        return "R"

    elif patient in NON_RESPONDERS:
        return "NR"

    return "Unknown"


def get_tissue_source(sample):

    if sample.endswith("_b"):
        return "Blood"

    elif sample.endswith("_t"):
        return "Tumor"

    return "Unknown"


def add_metadata(adata):

    adata.obs["Cancer type"] = "Endometrial Carcinoma"
    adata.obs["Tissue"] = "PBMC"

    adata.obs["Patient ID"] = (
        adata.obs_names
        .to_series()
        .apply(get_patient_id)
    )

    adata.obs["Response"] = (
        adata.obs_names
        .to_series()
        .apply(get_response)
    )


    adata.obs["Timepoint"] = (
        adata.obs_names
        .to_series()
        .apply(get_timepoint)
    )

    adata.obs["Sample ID"] = (
        adata.obs_names
        .to_series()
        .apply(get_sample_id)
    )

    adata.obs["Biopsy Site"] = (
        adata.obs_names
        .to_series()
        .apply(get_tissue_source)
    )

    adata.obs["Treatment"] = "Pembrolizumab"
    adata.obs["Treatment Class"] = "Anti-PD-1"

    print("-" * 80)

    unique_samples = (
        adata.obs[
            ["Sample ID", "Patient ID", "Response"]
        ]
        .drop_duplicates()
        .sort_values("Sample ID")
    )

    for _, row in unique_samples.iterrows():
        print(
            f"{row['Sample ID']:<15} "
            f"{row['Patient ID']:<8} "
            f"{row['Response']}"
        )

    print("-" * 80)

    return adata