
import subprocess
from pathlib import Path
import scanpy as sc


def load_data_rds(data_dir):
    data_dir = Path(data_dir).resolve()

    if not data_dir.is_dir():
        raise NotADirectoryError(
            f"RDS directory not found: {data_dir}"
        )

    r_script = Path(__file__).with_name("read_rds_and_convert_to_h5ad.R")

    rds_files = sorted(data_dir.glob("*.rds"))

    if not rds_files:
        raise FileNotFoundError(
            f"No .rds files found in {data_dir}"
        )

    datasets = []

    for file_path in rds_files:
        output_file = file_path.with_suffix(".h5ad")

        print(f"Loading {file_path.name}")

        subprocess.run(
            [
                "Rscript",
                str(r_script),
                str(file_path),
                str(output_file),
            ],
            check=True,
        )

        adata = sc.read_h5ad(output_file)
        datasets.append(adata)

    return datasets
