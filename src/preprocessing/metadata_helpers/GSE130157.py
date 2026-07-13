import pandas as pd

from ..config import DATASETS_TO_PREPROCESS




def add_metadata(adata):

    annotations = pd.read_csv(
        DATASETS_TO_PREPROCESS / "GSE130157/GSE130157.cell_annotations.txt",
        sep="\t"
    )

    annotations = annotations.set_index("Cell ID")

    common_cells = (
        adata.obs_names.intersection(
            annotations.index
        )
    )

    adata = adata[common_cells].copy()

    annotations = annotations.loc[common_cells]

    adata.obs["Patient ID"] = (
        annotations["Patient ID"].values
    )

    adata.obs["Timepoint"] = (
        annotations["Time Point"].values
    )

    adata.obs["Cancer type"] = (
        annotations["Cancer"].values
    )

    adata.obs["Treatment"] = (
        "mFOLFOX6 + Pembrolizumab"
    )

    adata.obs["Treatment Class"] = (
        "Anti-PD-1 + Chemotherapy"
    )

    adata.obs["Sample ID"] = (
        annotations["Patient ID"].astype(str)
        + "_"
        + annotations["Time Point"].astype(str)
    )

    adata.obs["Response"] = (
        annotations["Responder"]
        .replace(
            {
                "Responder": "R",
                "Non.Responder": "NR"
            }
        )
        .values
    )

    return adata