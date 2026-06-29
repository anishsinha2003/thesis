library(Seurat)
# read rds file (such as Liu et al.)
# remotes::install_github('cellgeni/sceasy')
# remotes::install_github('mojaveazure/seurat-disk')
find_h5ad_file <- function(folder) {
  # read input folder and find the list of name 
  file_name <- list.files(path = folder, pattern = "\\.rds$", full.names = TRUE)
  return(file_name)
}

convert_rds_to_h5ad <- function(list_files, out_folder) {
  for (file in list_files) {
    print(file)
    print("Read RDS file")
    seurat_obj <- readRDS(file)
    seurat_obj$dataset_path = file
    seurat_obj$dataset_alias = basename(file)
    # get sample id and sample name
    print(seurat_obj)
    # we need raw counts when transferring
    raw_counts <- GetAssayData(object = seurat_obj, assay = "RNA", layer = "counts")
    seurat_obj_raw <- CreateSeuratObject(counts = raw_counts, data = raw_counts)
    DefaultAssay(object = seurat_obj_raw) <- "RNA"
    sceasy::convertFormat(seurat_obj_raw, 
                          from="seurat", 
                          to="anndata", 
                          outFile=paste0(out_folder, "_", basename(file), ".h5ad"))
  }
}

# I need to do Liu et al. separately 
file_liu_et_al <- find_h5ad_file("~/Desktop/Benchmark-2025-Sep/external_validation_datasets/ji_et_al/liu_hcc/scRNAseqData_of_HCC")

convert_rds_to_h5ad(file_liu_et_al, "~/Desktop/ICI_Foundation_2026/ST_ICI_data/preprocessed/")

seurat_obj <- readRDS(file_liu_et_al)
# Cell type labelling
seurat_obj$Fib$cell_type <- "Fibroblasts"
seurat_obj$B$cell_type <- "B cells"
seurat_obj$Myeloid$cell_type <- "Myeloid"
seurat_obj$Endo$cell_type <- "Endothelial cells"
seurat_obj$T_NK$cell_type <- "T/NK cells"
seurat_obj$Hepato$cell_type <- "Hepatocytes"

# merge the seurat object
merged_obj_liu_et_al <- merge(
  x = seurat_obj[[1]], 
  y = seurat_obj[-1])

raw_counts <- GetAssayData(object = merged_obj_liu_et_al, assay = "RNA", layer = "counts")
merged_obj_liu_et_al_raw <- CreateSeuratObject(counts = raw_counts, data = raw_counts)




sceasy::convertFormat(merged_obj_liu_et_al_raw, 
                      from="seurat",       
                      to="anndata", 
                      outFile=paste0("merged_obj_liu_et_al.h5ad"))
###### Liu et al., GSE17994

file_gse179994 <- find_h5ad_file("~/Desktop/ICI_Foundation_2026/ST_ICI_data/GSE179994")

convert_rds_to_h5ad(file_gse179994, "~/Desktop/ICI_Foundation_2026/ST_ICI_data/preprocessed/")
file_gse179994_mtx <- readRDS(file_gse179994)
file_gse179994_obj <- CreateSeuratObject(counts = file_gse179994_mtx, data = file_gse179994_mtx)
DefaultAssay(object = file_gse179994_obj) <- "RNA"
sceasy::convertFormat(file_gse179994_obj, 
                      from="seurat",       
                      to="anndata", 
                      outFile=paste0("gse179994_obj.h5ad"))
