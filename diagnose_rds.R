
library(Seurat)
library(SeuratObject)

input_file <- commandArgs(trailingOnly = TRUE)[1]

obj <- readRDS(input_file)

cat("Loaded:", class(obj), "\n")

# Add the missing images slot if necessary
if (!"images" %in% slotNames(obj)) {
  slot(obj, "images") <- list()
}

cat("Checking object validity...\n")
print(tryCatch(
  validObject(obj),
  error = function(e) conditionMessage(e)
))

# Try updating the object to the current Seurat format
cat("Updating Seurat object...\n")
obj <- tryCatch(
  UpdateSeuratObject(obj),
  error = function(e) {
    message("Update failed: ", conditionMessage(e))
    obj
  }
)

cat("Checking RNA counts...\n")
counts <- tryCatch(
  GetAssayData(obj, assay = "RNA", layer = "counts"),
  error = function(e) {
    message("Could not extract counts: ", conditionMessage(e))
    NULL
  }
)

if (!is.null(counts)) {
  print(dim(counts))
  cat("RNA counts extracted successfully.\n")
}

# Save a copy only if the object can be validated
if (isTRUE(tryCatch(validObject(obj), error = function(e) FALSE))) {
  saveRDS(obj, sub("\\.rds$", "_updated.rds", input_file))
  cat("Saved updated object.\n")
} else {
  cat("Object is still invalid. No updated file saved.\n")
}
