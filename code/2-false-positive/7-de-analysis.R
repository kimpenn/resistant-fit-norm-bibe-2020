# Title     : DE Analysis for each downsampled group
# Created by: derek
# Created on: 7/9/20
library('Seurat')
library('ggplot2')
library('Matrix')
library('SingleCellExperiment')
library('scran')
library('scater')

PRJ_DIR <- "/home/derek/research/Kim-Lab/normalization-simulation/"
IN_DIR <- paste0(PRJ_DIR, "data/pbmc-33k/filtered_gene_bc_matrices/hg19/")
OUT_DIR <- paste0(PRJ_DIR, "exp/exp-12/out/")

file_list <- c("CD4_Memory", "CD4_Naive", "Mono_CD14", "B_Pre", "CD8_Memory",
               "CD8_Naive", "NK_Dim", "Mono_CD16", "CD8_Effector", "B_Pro",
               "DC", "NK_Bright", "Mk", "pDC")
WRT_FLG <- TRUE
# -----------------------------
# Export idents for scVI
# -----------------------------
for (i in seq(length(file_list))) {
  file.prefix <- file_list[i]
  file_name <- paste0(OUT_DIR, file.prefix, '-ident.Rds')
  ident <- readRDS(file_name)
  
  file_name <- paste0(OUT_DIR, file.prefix, '-ident.csv')
  write.csv(ident, file_name, row.names = FALSE)
}
# -----------------------------
# Resist Fit Normalization
# -----------------------------
num_marker <- c()
i <- 1

file.prefix <- file_list[11]
file_name <- paste0(OUT_DIR, file.prefix, 'cm-norm.mm.mtx')
print(file_name)
obs_resist_norm <- readMM(file_name)

file_name <- paste0(OUT_DIR, file.prefix, '-ident.Rds')
ident <- readRDS(file_name)

rownames(obs_resist_norm) <- seq(dim(obs_resist_norm)[1])
colnames(obs_resist_norm) <- seq(dim(obs_resist_norm)[2])
s <- CreateSeuratObject(counts = log1p(obs_resist_norm), min.cells = 5)
Idents(s) <- ident
markers <- FindMarkers(s, ident.1 = 'group1',
                       min.pct = 0.25,
                       slot = "data",
                       logfc.threshold = 0.3)

cat(dim(markers)[1], '\n')
# 3, 0, 5, 0, 2, 1, 4, 6, 2, 2, 22, 13, 34, 48  # logfc = 0.25
# 0, 0, 0, 0, 0, 0, 3, 2, 0, 0, 6, 4, 16, 18  # logfc = 0.3
# -----------------------------
# Log Normalization
# -----------------------------
num_marker <- c()
i <- 1

file.prefix <- file_list[14]
file_name <- paste0(OUT_DIR, file.prefix, '-cm.mm')
print(file_name)
obs_resist_norm <- readMM(file_name)

file_name <- paste0(OUT_DIR, file.prefix, '-ident.Rds')
ident <- readRDS(file_name)

rownames(obs_resist_norm) <- seq(dim(obs_resist_norm)[1])
colnames(obs_resist_norm) <- seq(dim(obs_resist_norm)[2])
s <- CreateSeuratObject(counts = obs_resist_norm, min.cells = 5)
s <- NormalizeData(s)
Idents(s) <- ident
markers <- FindMarkers(s, ident.1 = 'group1',
                       min.pct = 0.25,
                       slot = "data",
                       logfc.threshold = 0.3)

cat('[', i - 1, ']', ':', dim(markers)[1], '\n')
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 102, 92, 166, 305  # logfc = 0.25
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 57, 59, 119, 203  # logfc = 0.30
# -----------------------------
# scran
# -----------------------------
num_marker <- c()
i <- 1

file.prefix <- file_list[14]
#i <- i + 1
file_name <- paste0(OUT_DIR, file.prefix, '-cm.mm')
print(file_name)
obs_resist_norm <- readMM(file_name)

file_name <- paste0(OUT_DIR, file.prefix, '-ident.Rds')
ident <- readRDS(file_name)

rownames(obs_resist_norm) <- seq(dim(obs_resist_norm)[1])
colnames(obs_resist_norm) <- seq(dim(obs_resist_norm)[2])

s <- CreateSeuratObject(counts = obs_resist_norm, min.cells = 5)
sce <- SingleCellExperiment(assays = list(counts = s@assays$RNA@data))
#clusters <- quickCluster(sce)
sce <- computeSumFactors(sce)
sce <- logNormCounts(sce)
s <- CreateSeuratObject(counts = logcounts(sce))

Idents(s) <- ident
markers <- FindMarkers(s, ident.1 = 'group1',
                       min.pct = 0.25,
                       slot = "data",
                       logfc.threshold = 0.3)

cat(dim(markers)[1], '\n')
# 0, 0, 0, 0, 0, 0, 0,1, 0,1, 64, 26, 117, 197     # logfc = 0.25
# 0, 0, 0, 0, 0, 0, 0,1, 0,1, 27, 14, 74, 116     # logfc = 0.30
# -----------------------------
# SCTransform: corrected counts
# -----------------------------
file.prefix <- file_list[14]
file_name <- paste0(OUT_DIR, file.prefix, '-cm.mm')
print(file_name)
obs_resist_norm <- readMM(file_name)

file_name <- paste0(OUT_DIR, file.prefix, '-ident.Rds')
ident <- readRDS(file_name)

rownames(obs_resist_norm) <- seq(dim(obs_resist_norm)[1])
colnames(obs_resist_norm) <- seq(dim(obs_resist_norm)[2])
s <- CreateSeuratObject(counts = obs_resist_norm, min.cells = 5)
s <- SCTransform(s, do.scale = TRUE, do.center = TRUE,
                 return.only.var.genes = FALSE, clip.range = c(-Inf, Inf),
                 variable.features.n = Inf)
Idents(s) <- ident
markers <- FindMarkers(s, ident.1 = 'group1',
                       min.pct = 0.25,
                       slot = "data",
                       logfc.threshold = 0.3)

cat(dim(markers)[1], '\n')
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 9, 30  # logfc = 0.3
# -----------------------------
# SCTransform: Pearson Residual
# -----------------------------
file.prefix <- file_list[14]
file_name <- paste0(OUT_DIR, file.prefix, '-cm.mm')
print(file_name)
obs_resist_norm <- readMM(file_name)

file_name <- paste0(OUT_DIR, file.prefix, '-ident.Rds')
ident <- readRDS(file_name)

rownames(obs_resist_norm) <- seq(dim(obs_resist_norm)[1])
colnames(obs_resist_norm) <- seq(dim(obs_resist_norm)[2])
s <- CreateSeuratObject(counts = obs_resist_norm, min.cells = 5)
s <- SCTransform(s, do.scale = FALSE, do.center = FALSE,
                 return.only.var.genes = FALSE, clip.range = c(-Inf, Inf),
                 variable.features.n = Inf)
Idents(s) <- ident
markers <- FindMarkers(s, ident.1 = 'group1',
                       min.pct = 0.25,
                       slot = "scale.data",
                       logfc.threshold = 0.3)

cat(dim(markers)[1], '\n')
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 9, 1270  # logfc = 0.3

for (i in seq(14)) {
  file.prefix <- file_list[i]
  file_name <- paste0(OUT_DIR, file.prefix, '-cm.mm')
  print(file_name)
  obs_resist_norm <- readMM(file_name)
  print(dim(obs_resist_norm))
}