library(Matrix)

counts <- readRDS("/Users/anishsinha/Desktop/thesis/preprocessing/datasets_to_preprocess/original_raw_datasets/EGAS00001004809/1867-counts_cells_cohort2.rds")

writeMM(
  counts,
  "/Users/anishsinha/Desktop/thesis/preprocessing/datasets_to_preprocess/original_raw_datasets/EGAS00001004809/cohort2.mtx"
)

write.table(
  rownames(counts),
  "/Users/anishsinha/Desktop/thesis/preprocessing/datasets_to_preprocess/original_raw_datasets/EGAS00001004809/cohort2_genes.tsv",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  colnames(counts),
  "/Users/anishsinha/Desktop/thesis/preprocessing/datasets_to_preprocess/original_raw_datasets/EGAS00001004809/cohort2_barcodes.tsv",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)