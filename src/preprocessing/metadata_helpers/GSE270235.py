# ============================================================
# GSE270235 Metadata
# Metastatic Triple Negative Breast Cancer (mTNBC)
# Anti-PD-1 / PD-L1 + Chemotherapy
# ============================================================

import pandas as pd


# ------------------------------------------------------------
# Response
# ------------------------------------------------------------

def get_response(response):

    response = str(response)

    if response.lower() == "response":
        return "R"

    elif response.lower() == "progression":
        return "NR"

    return "Unknown"


# ------------------------------------------------------------
# Metadata
# ------------------------------------------------------------

def add_metadata(adata):

    adata.obs["Cancer type"] = (
        "Metastatic Triple Negative Breast Cancer"
    )

    adata.obs["Cancer Subtype"] = "TNBC"

    adata.obs["Disease Stage"] = "Metastatic"

    adata.obs["Tissue"] = "PBMC"

    adata.obs["Treatment"] = (
        "Anti-PD-1/PD-L1 + Chemotherapy"
    )

    adata.obs["Treatment Class"] = (
        "Immune Checkpoint Inhibitor + Chemotherapy"
    )

    # --------------------------------------------------------
    # Existing metadata from Seurat object
    # --------------------------------------------------------

    if "orig.ident" in adata.obs.columns:

        adata.obs["Sample ID"] = (
            adata.obs["orig.ident"]
            .astype(str)
        )

    else:

        adata.obs["Sample ID"] = "Unknown"

    if "patient" in adata.obs.columns:

        adata.obs["Patient ID"] = (
            "P"
            + adata.obs["patient"]
            .astype(str)
        )

    else:

        adata.obs["Patient ID"] = "Unknown"

    if "response" in adata.obs.columns:

        adata.obs["Response"] = (
            adata.obs["response"]
            .apply(get_response)
        )

    else:

        adata.obs["Response"] = "Unknown"

    if "timepoint" in adata.obs.columns:

        adata.obs["Timepoint"] = (
            adata.obs["timepoint"]
            .astype(str)
        )

    else:

        adata.obs["Timepoint"] = "Unknown"

    return adata