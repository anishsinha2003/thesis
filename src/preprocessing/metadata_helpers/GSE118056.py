
import pandas as pd


def add_metadata(adata):
    barcodes = pd.Index(adata.obs_names.astype(str))

    # GEO documents -1 as PBMC and -2 as tumour.
    suffix = barcodes.str.extract(r"-(\d+)(?:-\d+)?$", expand=False)

    tissue_map = {
        "1": "PBMC",
        "2": "Tumor",
    }

    adata.obs["Tissue"] = suffix.map(tissue_map)

    adata.obs["Cancer_Type"] = "Merkel Cell Carcinoma"
    adata.obs["Time_Point"] = adata.obs["Tissue"].map({
        "PBMC": "Relapse_D559",
        "Tumor": "Relapse_D565",
    })

    # Both samples are from a relapse/resistance time point.
    # Do not treat these as separate responder/non-responder labels.
    adata.obs["Response"] = "Relapse"

    if adata.obs["Tissue"].isna().any():
        print(
            "Warning: Some barcodes could not be assigned a tissue. "
            "Check the barcode suffix format."
        )

    return adata
