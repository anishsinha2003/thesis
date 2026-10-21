
library(Seurat)
library(sceasy)

# Find RDS files in the dataset folder
find_h5ad_file <- function(folder) {
  file_name <- list.files(
    path = folder,
    pattern = "\\.rds$",
    full.names = TRUE
  )
  return(file_name)
}

# Convert RDS files to H5AD
convert_rds_to_h5ad <- function(list_files, out_folder) {
  for (file in list_files) {
    print(file)
    print("Read RDS file")

    seurat_obj <- readRDS(file)

    seurat_obj$dataset_path <- file
    seurat_obj$dataset_alias <- basename(file)

    print(seurat_obj)

    # Extract raw RNA counts
    raw_counts <- GetAssayData(
      object = seurat_obj,
      assay = "RNA",
      layer = "counts"
    )

    # Create a Seurat object containing raw counts
    seurat_obj_raw <- CreateSeuratObject(
      counts = raw_counts,
      data = raw_counts
    )

    DefaultAssay(object = seurat_obj_raw) <- "RNA"

    # Convert to H5AD
    sceasy::convertFormat(
      seurat_obj_raw,
      from = "seurat",
      to = "anndata",
      outFile = file.path(
        out_folder,
        paste0(tools::file_path_sans_ext(basename(file)), ".h5ad")
      )
    )
  }
}

# ============================================================
# GSE179784
# ============================================================

input_folder <- file.path(
  "/Users/anishsinha/Desktop/thesis/preprocessing",
  "datasets_to_preprocess",
  "GSE179784"
)

output_folder <- input_folder

# Find the RDS file
file_gse179784 <- find_h5ad_file(input_folder)

if (length(file_gse179784) == 0) {
  stop("No RDS files found in: ", input_folder)
}

# Create H5AD file(s)
convert_rds_to_h5ad(
  file_gse179784,
  output_folder
)
