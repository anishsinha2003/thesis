
# ============================================================
# Convert GSE179784 RDS file to H5AD
# ============================================================

# Install missing packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}

for (pkg in c("SingleCellExperiment", "zellkonverter")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)

# Input folder containing the RDS file
input_dir <- "/Users/anishsinha/Desktop/thesis/preprocessing/datasets_to_preprocess/GSE179784"

# Find RDS files in the folder
rds_files <- list.files(
  input_dir,
  pattern = "\\.rds$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(rds_files) == 0) {
  stop("No .rds files found in: ", input_dir)
}

if (length(rds_files) > 1) {
  stop(
    "Multiple RDS files found. Please choose the correct one:\n",
    paste(basename(rds_files), collapse = "\n")
  )
}

rds_file <- rds_files[1]

# Output H5AD file in the same folder
output_file <- file.path(input_dir, "GSE179784.h5ad")

cat("Reading:", rds_file, "\n")

# Load the RDS object
obj <- readRDS(rds_file)
obj <- SeuratObject::UpdateSeuratObject(obj)

cat("Object class:", paste(class(obj), collapse = ", "), "\n")

# Convert supported object types
if (inherits(obj, "Seurat")) {

  # Prefer the RNA assay if available
  if ("RNA" %in% SeuratObject::Assays(obj)) {
    SeuratObject::DefaultAssay(obj) <- "RNA"
  }

  sce <- Seurat::as.SingleCellExperiment(obj)

} else if (inherits(obj, "SingleCellExperiment")) {

  sce <- obj

} else if (inherits(obj, "SummarizedExperiment")) {

  sce <- SingleCellExperiment::SingleCellExperiment(obj)

} else if (is.matrix(obj) || inherits(obj, "Matrix")) {

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = obj)
  )

} else if (is.data.frame(obj)) {

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = as.matrix(obj))
  )

} else {

  stop(
    "Unsupported RDS object type: ",
    paste(class(obj), collapse = ", "),
    ". Inspect the object before converting."
  )
}

# Display dimensions
cat("Genes:", nrow(sce), "\n")
cat("Cells:", ncol(sce), "\n")

if (nrow(sce) == 0 || ncol(sce) == 0) {
  stop("The object has zero genes or zero cells.")
}

# Export to H5AD
zellkonverter::writeH5AD(
  sce,
  file = output_file
)

cat("\nConversion complete!\n")
cat("Saved to:", output_file, "\n")
