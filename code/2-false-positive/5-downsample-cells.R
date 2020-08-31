# Title     : Split and Downsample each cluster
# Created by: derek
# Created on: 7/8/20
library('Seurat')
library('Matrix')
library('dplyr')
library('parallel')

PRJ_DIR <- "/home/derek/research/Kim-Lab/normalization-simulation/"
IN_DIR <- paste0(PRJ_DIR, "data/pbmc-33k/filtered_gene_bc_matrices/hg19/")
OUT_DIR <- paste0(PRJ_DIR, "exp/exp-12/out/")

WRT_FLG <- TRUE
rand_seed <- 42 # Don't Panic!
downsample.rate <- 0.5
# -----------------------------
# Read data and Preprocess
# -----------------------------
# load the 33k PBMC data
cm <- Read10X(data.dir = IN_DIR)
rownames(cm) <- gsub('_', '-', make.unique(make.names(rownames(cm))))
# filter out very highly or lowly sequenced cells
log.umi <- log10(colSums(cm))
keep <- abs(scale(log.umi)[, 1]) <= 2
cm <- cm[, keep]
# load cell identities (based on a previous clustering)
file_name <- paste0(IN_DIR, 'pbmc33k_identities.Rds')
ident <- readRDS(file_name)
cell_type_list <- as.character(unique(ident))
summary(ident)
# [1] "CD4_Memory"   "CD4_Naive"    "Mono_CD14"    "B_Pre"        "CD8_Memory"
# [6] "CD8_Naive"    "NK_Dim"       "Mono_CD16"    "CD8_Effector" "B_Pro"
#[11] "DC"           "NK_Bright"    "Mk"           "pDC"

# -----------------------------
# Downsample and crate the count matrix
# -----------------------------
for (i in seq(length(cell_type_list))) {
#for (i in seq(1)) {
  SELECTED_CELL_TYPE <- cell_type_list[i]
  cat('working on ', SELECTED_CELL_TYPE, '...\n')
  # select barcodes by the cell type
  sel1 <- names(ident)[ident == SELECTED_CELL_TYPE]
  sel1 <- sel1[sel1 %in% colnames(cm)]
  sel2 <- names(ident)[ident == SELECTED_CELL_TYPE]
  sel2 <- sel2[sel2 %in% colnames(cm)]
  # divide the cells up into two groups
  set.seed(i)
  sel1 <- sample(sel1, ceiling(length(sel1) / 2))
  sel2 <- setdiff(sel2, sel1)
  cm1 <- cm[, sel1]
  cm2 <- cm[, sel2]
  cat('dim(cm1):', dim(cm1), '\n')
  cat('dim(cm2):', dim(cm2), '\n')
  cat('\t Downsampling...\n')
  cm2 <- SampleUMI(cm2, max.umi = floor(colSums(cm2) * downsample.rate))
  colnames(cm2) <- paste0(colnames(cm2), '_ds')
  # combine the two matrices
  this.cm <- cbind(cm1, cm2)
  this.idents <- factor(c(rep('group1', ncol(cm1)), rep('group2', ncol(cm2))),
                        levels = c('group1', 'group2'), ordered = TRUE)

  if (WRT_FLG) {
    file_name <- paste0(OUT_DIR, SELECTED_CELL_TYPE, '-cm.mm')
    if (!file.exists(file_name)) {
      cat('Writing to: ', file_name, '\n')
      #write.csv(file = file_name, as.matrix(cm))
      writeMM(this.cm, file = file_name)
    }

    file_name <- paste0(OUT_DIR, SELECTED_CELL_TYPE, '-ident.Rds')
    if (!file.exists(file_name)) {
      cat('Writing to: ', file_name, '\n')
      saveRDS(this.idents, file = file_name)
    }
  }
}
