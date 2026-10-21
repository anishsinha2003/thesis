import pandas as pd


def add_metadata(adata):
    # Cell Ranger aggregated matrices may use suffixes
    # to distinguish the original libraries.
    barcodes = pd.Index(adata.obs_names.astype(str))

    sample_suffix = barcodes.str.extract(r"-(\d+)$", expand=False)

    sample_map = {
        "1": "Tumor_PreRx",
        "2": "Tumor_Acquired_Resistance",
    }

    adata.obs["Sample_ID"] = sample_suffix.map(sample_map)

    adata.obs["Patient_ID"] = "2586-4"
    adata.obs["Tissue"] = "Tumor"
    adata.obs["Cancer_Type"] = "Merkel Cell Carcinoma"

    adata.obs["Time_Point"] = adata.obs["Sample_ID"].map({
        "Tumor_PreRx": "PreRx",
        "Tumor_Acquired_Resistance": "Day+615",
    })

    # These are longitudinal time points, not responder/non-responder
    # patient labels. Do not treat them as equivalent to R/NR.
    adata.obs["Response"] = "Unknown"

    # Warn if the assumed barcode mapping did not work.
    if adata.obs["Sample_ID"].isna().any():
        print(
            "Warning: Some cell barcodes were not mapped. "
            "Inspect adata.obs_names and confirm the sample labels."
        )

    return adata