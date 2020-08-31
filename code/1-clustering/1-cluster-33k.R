# Title     : Create Monocytes CD 14 data from 33k pbmc dataset
# Created by: derek
# Created on: 7/2/20
# Reference: SCTransform: https://osf.io/kc9dg/
library('Seurat')
library('ggplot2')
library('Matrix')
library('dplyr')
library('parallel')
library('loomR')
library('patchwork')
library('SingleCellExperiment')
library('scran')
library('scater')
PRJ_DIR <- "/home/derek/research/Kim-Lab/normalization-simulation/"
IN_DIR <- paste0(PRJ_DIR, "data/pbmc-33k/filtered_gene_bc_matrices/hg19/")
OUT_DIR <- paste0(PRJ_DIR, "exp/exp-12/out/")

WRT_FLG <- TRUE
SELECTED_CELL_TYPE <- 'Mono_CD14'
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
cm <- cm[, keep]      # 32738 31080
# load cell identities (based on a previous clustering)
file_name <- paste0(IN_DIR, 'pbmc33k_identities.Rds')
celltype <- readRDS(file_name)
# -----------------------------      
# Export count matrix for resistant fit
# -----------------------------
if (WRT_FLG) {
  file_name <- paste0(OUT_DIR, 'cm.csv')
  if (!file.exists(file_name)) {
    cat('Writing to: ', file_name, '\n')
    #write.csv(file = file_name, as.matrix(cm))
    writeMM(cm, file = file_name)
  }
}
# -----------------------------
# Export count matrix for scVI
# -----------------------------
if (WRT_FLG) {
  file_name <- paste0(OUT_DIR, 'pbmc33k_identities.csv')
  write.csv(celltype, file_name)
}

# -----------------------------
# Clustering without Normalization
# -----------------------------
s <- CreateSeuratObject(counts = cm, min.cells = 5)
s <- FindVariableFeatures(s, selection.method = "vst", nfeatures = 2000)
s <- ScaleData(s)
s[['RNA']]@scale.data[s[['RNA']]@scale.data < -10] <- -10
s[['RNA']]@scale.data[s[['RNA']]@scale.data > 10] <- 10
s <- RunPCA(s, verbose = FALSE, npcs = 25, features = VariableFeatures(object = s))
#s <- RunPCA(s, verbose = FALSE, npcs = 20, features = rownames(s[['RNA']]@scale.data))
s <- RunUMAP(s, dims = 1:25)
s <- FindNeighbors(object = s, dims = 1:25, verbose = FALSE)
s <- FindClusters(object = s, resolution = 0.6, verbose = FALSE)
DimPlot(s, reduction = "umap")
s@meta.data$celltype <- NA
s@meta.data[names(celltype), 'celltype'] <- as.character(celltype)
tmp1 <- table(as.character(s@meta.data$seurat_clusters), as.character(s@meta.data$celltype))
print(tmp1)
# -----------------------------
# Clustering with Log Norm
# -----------------------------
s <- CreateSeuratObject(counts = cm, min.cells = 5)
s <- NormalizeData(s)
s <- FindVariableFeatures(s, selection.method = "vst", nfeatures = 2000)
s <- ScaleData(s)
s[['RNA']]@scale.data[s[['RNA']]@scale.data < -10] <- -10
s[['RNA']]@scale.data[s[['RNA']]@scale.data > 10] <- 10
s <- RunPCA(s, verbose = FALSE, npcs = 25, features = VariableFeatures(object = s))
#s <- RunPCA(s, verbose = FALSE, npcs = 20, features = rownames(s[['RNA']]@scale.data))
s <- RunUMAP(s, dims = 1:25)
s <- FindNeighbors(object = s, dims = 1:25, verbose = FALSE)
s <- FindClusters(object = s, resolution = 0.6, verbose = FALSE)
DimPlot(s, reduction = "umap")
ident.cluster <- Idents(s)

s@meta.data$celltype <- NA
s@meta.data[names(celltype), 'celltype'] <- as.character(celltype)
tmp1 <- table(as.character(s@meta.data$seurat_clusters), as.character(s@meta.data$celltype))
print(tmp1)

if (WRT_FLG) {
  file_name <- paste0(OUT_DIR, 'logNorm.Rds')
  if (!file.exists(file_name)) {
    cat('Writing to: ', file_name, '\n')
    res <- list(md = s@meta.data, umap = s[['umap']])
    saveRDS(res, file = file_name)
  }
}
# -----------------------------
# Clustering with SCTransform
# -----------------------------
s <- CreateSeuratObject(counts = cm, min.cells = 5)
s <- SCTransform(s, do.scale = TRUE, do.center = TRUE,
                 return.only.var.genes = FALSE, clip.range = c(-Inf, Inf),
                 variable.features.n = Inf)
s[['RNA']]@scale.data[s[['RNA']]@scale.data < -10] <- -10
s[['RNA']]@scale.data[s[['RNA']]@scale.data > 10] <- 10
s <- RunPCA(s, verbose = FALSE, npcs = 25, features = rownames(s[['RNA']]@scale.data))
s <- RunUMAP(s, dims = 1:25)

s <- FindNeighbors(object = s, dims = 1:25, verbose = FALSE)
s <- FindClusters(object = s, resolution = 0.6, verbose = FALSE)
DimPlot(s, reduction = "umap")
ident.cluster <- Idents(s)

s@meta.data$celltype <- NA
s@meta.data[names(celltype), 'celltype'] <- as.character(celltype)
tmp1 <- table(as.character(s@meta.data$seurat_clusters), as.character(s@meta.data$celltype))
print(tmp1)

if (WRT_FLG) {
  file_name <- paste0(OUT_DIR, 'SCTransform.Rds')
  if (!file.exists(file_name)) {
    cat('Writing to: ', file_name, '\n')
    res <- list(md = s@meta.data, umap = s[['umap']])
    saveRDS(res, file = file_name)
  }
}
ident <- c('Monocytes', 'T cells', 'B cells', 'T cells', 'T cells', 'T cells',   # 0-5
           'Monocytes', 'NK cells', 'T cells', 'T cells', 'T cells',  # 6-10
           'NA', 'NA', 'NA', 'NA', 'NA', 'NA')   # 11-16
levels(ident.cluster) <- ident
Idents(s) <- ident.cluster
DimPlot(s, reduction = "umap", cells = colnames(s)[Idents(s) != 'NA'])
# -----------------------------
# Clustering with Resistant Fit
# -----------------------------
#file_name <- paste0(OUT_DIR, 'norm.cm-iter-1-no-intercept.mm.mtx')
file_name <- paste0(OUT_DIR, 'norm.cm-iter-1-intercept.mm.mtx')
obs_resist_norm <- readMM(file_name)
rownames(obs_resist_norm) <- rownames(cm)
colnames(obs_resist_norm) <- colnames(cm)
obs_resist_norm[obs_resist_norm < 0] <- 0

s <- CreateSeuratObject(counts = log1p(10 * obs_resist_norm), min.cells = 5)
s <- FindVariableFeatures(s, selection.method = "vst", nfeatures = 2000)
s <- ScaleData(s)
s[['RNA']]@scale.data[s[['RNA']]@scale.data < -10] <- -10
s[['RNA']]@scale.data[s[['RNA']]@scale.data > 10] <- 10
s <- RunPCA(s, verbose = FALSE, npcs = 25, features = VariableFeatures(object = s))
#s <- RunPCA(s, verbose = FALSE, npcs = 20, features = rownames(s[['RNA']]@scale.data))
s <- RunUMAP(s, dims = 1:25)
s <- FindNeighbors(object = s, dims = 1:25, verbose = FALSE)
s <- FindClusters(object = s, resolution = 0.6, verbose = FALSE)
DimPlot(s, reduction = "umap")

ident.cluster <- Idents(s)
ident.cluster.bk <- Idents(s)
s@meta.data$celltype <- NA
s@meta.data[names(celltype), 'celltype'] <- as.character(celltype)
tmp1 <- table(as.character(s@meta.data$seurat_clusters), as.character(s@meta.data$celltype))
print(tmp1)

if (WRT_FLG) {
  file_name <- paste0(OUT_DIR, 'ResistantFit.Rds')
  if (!file.exists(file_name)) {
    cat('Writing to: ', file_name, '\n')
    res <- list(md = s@meta.data, umap = s[['umap']])
    saveRDS(res, file = file_name)
  }
}

ident <- c('T cells', 'Monocytes', 'T cells', 'B cells', 'T cells', 'Monocytes',   # 0-5
           'NK cells', 'T cells', 'T cells', 'T cells', 'T cells',   # 6-10
           'NA', 'Monocytes', 'B cells', 'NA', 'NA')   # 11-15
levels(ident.cluster) <- ident
Idents(s) <- ident.cluster
DimPlot(s, reduction = "umap", cells = colnames(s)[Idents(s) != 'NA'])
# -----------------------------
# Clustering with SCRAN
# -----------------------------
s <- CreateSeuratObject(counts = cm, min.cells = 5)
sce <- SingleCellExperiment(assays = list(counts = s@assays$RNA@data))
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters = clusters)
sce <- logNormCounts(sce)

s <- CreateSeuratObject(counts = logcounts(sce))
s <- FindVariableFeatures(s, selection.method = "vst", nfeatures = 2000)
s <- ScaleData(s)
s[['RNA']]@scale.data[s[['RNA']]@scale.data < -10] <- -10
s[['RNA']]@scale.data[s[['RNA']]@scale.data > 10] <- 10
s <- RunPCA(s, verbose = FALSE, npcs = 25, features = VariableFeatures(object = s))
#s <- RunPCA(s, verbose = FALSE, npcs = 20, features = rownames(s[['RNA']]@scale.data))
s <- RunUMAP(s, dims = 1:25)
s <- FindNeighbors(object = s, dims = 1:25, verbose = FALSE)
s <- FindClusters(object = s, resolution = 0.6, verbose = FALSE)
DimPlot(s, reduction = "umap")

ident.cluster <- Idents(s)
s@meta.data$celltype <- NA
s@meta.data[names(celltype), 'celltype'] <- as.character(celltype)
tmp1 <- table(as.character(s@meta.data$seurat_clusters), as.character(s@meta.data$celltype))
print(tmp1)

if (WRT_FLG) {
  file_name <- paste0(OUT_DIR, 'SCRAN.Rds')
  if (!file.exists(file_name)) {
    cat('Writing to: ', file_name, '\n')
    res <- list(md = s@meta.data, umap = s[['umap']])
    saveRDS(res, file = file_name)
  }
}

ident <- c('T cells', 'Monocytes', 'B cells', 'T cells', 'T cells', 'Monocytes',   # 0-5
           'NK cells', 'T cells', 'T cells', 'T cells', 'T cells', 'NA', 'Monocytes',
           'B cells', 'NA', 'NA', 'NA')   # 6-14
levels(ident.cluster) <- ident
Idents(s) <- ident.cluster
DimPlot(s, reduction = "umap", cells = colnames(s)[Idents(s) != 'NA'])

# -----------------------------
# Clustering with scVI
# -----------------------------
file_name <- paste0(OUT_DIR, 'scvi-embedding.csv')
scvi.df <- read.csv(file_name)
rownames(scvi.df) <- scvi.df$barcode
tmp1 <- table(as.character(scvi.df$leiden), as.character(celltype[scvi.df$barcode]))
print(tmp1)