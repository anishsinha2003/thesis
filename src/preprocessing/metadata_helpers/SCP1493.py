
import pandas as pd
from pathlib import Path
from ..config import DATASETS_TO_PREPROCESS

def add_metadata(adata):
    input_dir = Path(DATASETS_TO_PREPROCESS) / "SCP1493"

    metadata = pd.read_csv(
        input_dir / "BM_c1000_base_metadata.tsv",
        sep="\t",
    )

    if "NAME" not in metadata.columns:
        raise ValueError(
            "Metadata has no NAME column. "
            f"Available columns: {metadata.columns.tolist()}"
        )

    metadata["NAME"] = metadata["NAME"].astype(str)
    metadata = metadata.drop_duplicates("NAME").set_index("NAME")

    # Match metadata rows to the cells in the expression matrix.
    matched = metadata.reindex(adata.obs_names)

    for column in matched.columns:
        adata.obs[column] = matched[column].to_numpy()

    missing = matched.index[matched.isna().all(axis=1)]

    if len(missing) > 0:
        print(
            f"Warning: metadata was not matched for "
            f"{len(missing)} cells. Check the NAME identifiers."
        )

    return adata
