# Title     : subset CD4 memory T and CD8 memory T cell
# Created by: derek
# Created on: 7/10/20

library('Seurat')
library('SingleCellExperiment')
library('scran')
library('scater')
library('Matrix')

PRJ_DIR <- "/home/derek/research/Kim-Lab/normalization-simulation/"
IN_DIR <- paste0(PRJ_DIR, "data/pbmc-33k/filtered_gene_bc_matrices/hg19/")
OUT_DIR <- paste0(PRJ_DIR, "exp/exp-14/out/")

WRT_FLG <- TRUE
#SELECTED_CELL_TYPE <- c('CD4_Memory', 'CD8_Memory')
SELECTED_CELL_TYPE <- c('CD4_Memory', 'CD8_Effector')
rand_seed <- 42 # Don't Panic!
# -----------------------------
# Read data and Preprocess
# -----------------------------
# load the 33k PBMC data
cm <- Read10X(data.dir = IN_DIR)
rownames(cm) <- gsub('_', '-', make.unique(make.names(rownames(cm))))
# filter out very highly or lowly sequenced cells
log.umi <- log10(colSums(cm))
keep <- abs(scale(log.umi)[, 1]) <= 2
cm <- cm[, keep]      # 32738 31080
# load cell identities (based on a previous clustering)
file_name <- paste0(IN_DIR, 'pbmc33k_identities.Rds')
celltype <- readRDS(file_name)
# -----------------------------
# Subset Cells
# -----------------------------
s <- Seurat::CreateSeuratObject(counts = cm, min.cells = 5)
s <- s[, names(celltype)]
Idents(s) <- celltype[colnames(s)]
s.cd4.cd8 <- subset(s, idents = SELECTED_CELL_TYPE)
# -----------------------------
# Export data
# -----------------------------
if (WRT_FLG) {
  # Rds for SCRAN, Log Norm and SCTransform
  #file_name <- paste0(OUT_DIR, 's.cd4.cd8.Rds')
  file_name <- paste0(OUT_DIR, 's.cd4.cd8.eff.Rds')
  if (!file.exists(file_name)) {
    cat('Writing to: ', file_name, '\n')
    saveRDS(s.cd4.cd8, file = file_name)
  }
  # matrix market for RFnorm
  #file_name <- paste0(OUT_DIR, 's.cd4.cd8.mm.mtx')
  file_name <- paste0(OUT_DIR, 's.cd4.cd8.eff.mm.mtx')
  if (!file.exists(file_name)) {
    cat('Writing to: ', file_name, '\n')
    #write.csv(file = file_name, as.matrix(cm))
    writeMM(s.cd4.cd8@assays$RNA@data, file = file_name)
  }
  # csv and barcode for scanpy and scvi
  #file_name <- paste0(OUT_DIR, 's.cd4.cd8.csv')
  file_name <- paste0(OUT_DIR, 's.cd4.cd8.eff.csv')
  if (!file.exists(file_name)) {
    cat('Writing to: ', file_name, '\n')
    write.csv(s.cd4.cd8@assays$RNA@data, file_name, quote = FALSE)
  }
  #file_name <- paste0(OUT_DIR, 's.cd4.cd8.barcode.csv')
  file_name <- paste0(OUT_DIR, 's.cd4.cd8.eff.barcode.csv')
  if (!file.exists(file_name)) {
    cat('Writing to: ', file_name, '\n')
    write.csv(Idents(s.cd4.cd8), file_name, quote = FALSE)
  }
}