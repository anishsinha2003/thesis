# ============================================================
# GSE243572 Metadata
# Hepatocellular Carcinoma
# Regorafenib + Nivolumab
# ============================================================
import re
def get_response(sample):
    sample = str(sample)

    if "LPFS" in sample:
        return "R"
    elif "SPFS" in sample:
        return "NR"

    return "Unknown"


def get_timepoint(sample):
    sample = str(sample)

    if "C1D1" in sample:
        return "Pre"
    elif "C2D1" in sample:
        return "On Treatment"
    elif "C3D1" in sample:
        return "On Treatment"

    return "Unknown"


def extract_sample_id(sample):

    sample = str(sample)

    match = re.search(r"(LPFS|SPFS).*?(?=_GEX)", sample)

    if match:
        return match.group(0)

    return sample

def extract_patient_id(sample):

    sample = str(sample)

    match = re.search(
        r"(?:LPFS|SPFS)_[^_]+_(.*?)(?:_GEX|_TCR)",
        sample
    )

    if match:
        return match.group(1)

    return "Unknown"


def add_metadata(adata):

    adata.obs["Cancer type"] = "Hepatocellular Carcinoma"
    adata.obs["Tissue"] = "PBMC"

    adata.obs["Treatment"] = "Regorafenib + Nivolumab"
    adata.obs["Treatment Class"] = (
        "Anti-PD-1 + Tyrosine Kinase Inhibitor"
    )

    adata.obs["Response"] = (
        adata.obs["sample_id"]
        .apply(get_response)
    )

    adata.obs["Sample ID"] = (
        adata.obs["sample_id"]
        .apply(extract_sample_id)
    )

    adata.obs["Patient ID"] = (
        adata.obs["sample_id"]
        .apply(extract_patient_id)
    )

    adata.obs["Timepoint"] = (
        adata.obs["sample_id"]
        .apply(get_timepoint)
    )

    return adata