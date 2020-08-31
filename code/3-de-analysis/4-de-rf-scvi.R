# Title     : Get de genes for resistant fit and scvi
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
file_name <- paste0(OUT_DIR, 's.cd4.cd8.eff.Rds')
s.cd4.cd8 <- readRDS(file = file_name)
# -----------------------------
# DE: Resistant Fit Normalization
# -----------------------------
file_name <- paste0(OUT_DIR, 's.norm.cd4.cd8.eff.mm.mtx')
obs_resist_norm <- readMM(file_name)

rownames(obs_resist_norm) <- rownames(s.cd4.cd8)
colnames(obs_resist_norm) <- colnames(s.cd4.cd8)

norm.resist <- Seurat::CreateSeuratObject(counts = log1p(obs_resist_norm),
                                          min.cells = 5)
Idents(norm.resist) <- Idents(s.cd4.cd8)
markers <- FindMarkers(norm.resist,
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
  file_name <- paste0(OUT_DIR, '2-marker.rf.cd4.cd8.eff.csv')
  if (!file.exists(file_name)) {
    cat('Writing to: ', file_name, '\n')
    write.csv(markers, file = file_name)
  }
}
# -----------------------------
# DE: SCVI
# -----------------------------
file_name <- paste0(OUT_DIR, '2-marker.scvi.cd4.cd8.eff.full.csv')
markers <- read.csv(file = file_name, row.names=1)
markers['symbol'] <- row.names(markers)
select.idx <- abs(markers['bayes_factor']) > 3.0
if (WRT_FLG) {
  # Rds for SCRAN, Log Norm and SCTransform
  file_name <- paste0(OUT_DIR, '2-marker.scvi.cd4.cd8.eff.csv')
  if (!file.exists(file_name)) {
    cat('Writing to: ', file_name, '\n')
    write.csv(markers[select.idx,c('bayes_factor', 'symbol')], file = file_name)
  }
}
