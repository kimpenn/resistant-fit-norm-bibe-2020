# Title     : DE analysis for LogNorm, SCRAN and SCTransform
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
#file_name <- paste0(OUT_DIR, 's.cd4.cd8.Rds')
file_name <- paste0(OUT_DIR, 's.cd4.cd8.eff.Rds')
s.cd4.cd8 <- readRDS(file = file_name)
# -----------------------------
# Log Norm and DE
# -----------------------------
s.cd4.cd8.log <- subset(s.cd4.cd8)
markers <- FindMarkers(s.cd4.cd8.log,
                       ident.1 = SELECTED_CELL_TYPE[1],
                       min.pct = 0.25,
                       slot = "data",
                       logfc.threshold = 0.3)
markers['symbol'] <- rownames(markers)
n.de <- dim(markers)[1]
n.de.pos <- sum(markers$avg_logFC > 0)
n.de.neg <- sum(markers$avg_logFC < 0)
if (WRT_FLG) {
  # Rds for SCRAN, Log Norm and SCTransform
  file_name <- paste0(OUT_DIR, '2-marker.log.cd4.cd8.eff.csv')
  if (!file.exists(file_name)) {
    cat('Writing to: ', file_name, '\n')
    write.csv(markers, file = file_name)
  }
}

# -----------------------------
# SCTransform and DE
# -----------------------------
s.cd4.cd8.sct <- subset(s.cd4.cd8)
s.cd4.cd8.sct <- SCTransform(s.cd4.cd8.sct, do.scale = TRUE, do.center = TRUE,
                             return.only.var.genes = FALSE, clip.range = c(-Inf, Inf),
                             variable.features.n = Inf)
markers <- FindMarkers(s.cd4.cd8.sct,
                       ident.1 = SELECTED_CELL_TYPE[1],
                       min.pct = 0.25,
                       slot = "data",
                       logfc.threshold = 0.3)
markers['symbol'] <- rownames(markers)
n.de <- dim(markers)[1]
n.de.pos <- sum(markers$avg_logFC > 0)
n.de.neg <- sum(markers$avg_logFC < 0)
if (WRT_FLG) {
  # Rds for SCRAN, Log Norm and SCTransform
  file_name <- paste0(OUT_DIR, '2-marker.sct.cd4.cd8.eff.csv')
  if (!file.exists(file_name)) {
    cat('Writing to: ', file_name, '\n')
    write.csv(markers, file = file_name)
  }
}

# -----------------------------
# SCRAN and DE
# -----------------------------
sce <- SingleCellExperiment(assays = list(counts = s.cd4.cd8@assays$RNA@data))
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters = clusters)
sce <- logNormCounts(sce)
s.cd4.cd8.scran <- CreateSeuratObject(counts = logcounts(sce))
Idents(s.cd4.cd8.scran) <- Idents(s.cd4.cd8)
markers <- FindMarkers(s.cd4.cd8.scran,
                       ident.1 = SELECTED_CELL_TYPE[1],
                       min.pct = 0.25,
                       slot = "data",
                       logfc.threshold = 0.3)
markers['symbol'] <- rownames(markers)
n.de <- dim(markers)[1]
n.de.pos <- sum(markers$avg_logFC > 0)
n.de.neg <- sum(markers$avg_logFC < 0)
if (WRT_FLG) {
  # Rds for SCRAN, Log Norm and SCTransform
  file_name <- paste0(OUT_DIR, '2-marker.scran.cd4.cd8.eff.csv')
  if (!file.exists(file_name)) {
    cat('Writing to: ', file_name, '\n')
    write.csv(markers, file = file_name)
  }
}