# import scanpy as sc
# from pathlib import Path


# def combine_processed_data(processed_dir):


#     files = sorted(
#         Path(processed_dir).glob("*.h5ad")
#     )

#     print("Files Found:")
#     for file in files:
#         print(file)

#     adatas = []

#     for file in files:

#         adata = sc.read_h5ad(file)

#         # Make cell barcodes unique across samples
#         adata.obs_names = (
#             file.stem + "_" + adata.obs_names.astype(str)
#         )

#         print(
#             f"Loaded {file.name}: "
#             f"{adata.n_obs} cells x {adata.n_vars} genes"
#         )

#         adatas.append(adata)

#     combined = sc.concat(
#         adatas,
#         join="outer",
#         merge="same"
#     )

#     print("\nCombined Successfully")

#     return combined