# Title     : Visulaze the normalization results
# Created by: derek
# Created on: 7/7/20
library('Matrix')
library('ggplot2')
library('reshape2')
library('scales')
library('dplyr')
library('cowplot')
library(scater)

PRJ_DIR <- "/home/derek/research/Kim-Lab/normalization-simulation/"
IN_DIR <- paste0(PRJ_DIR, "data/pbmc-33k/filtered_gene_bc_matrices/hg19/")
OUT_DIR <- paste0(PRJ_DIR, "exp/exp-12/out/")
WRT_FLG <- TRUE
file_list <- c('logNorm.Rds', 'SCTransform.Rds', 'ResistantFit.Rds', 'SCRAN.Rds')

# -----------------------------
# Hardcod the celltypes of each clusters
# -----------------------------
celltype <- list(
  'log-norm' = c('T cells', 'Monocytes', 'T cells', 'B cells', 'T cells', 'Monocytes',   # 0-5
                 'T cells', 'NK cells', 'T cells', 'T cells', 'NA',  # 6 - 10
                 'NA', 'NA', 'NA', 'NA', 'NA'),   # 11 - 15
  'SCTransform' = c('Monocytes', 'T cells', 'B cells', 'T cells', 'T cells', 'T cells',   # 0-5
                    'Monocytes', 'NK cells', 'T cells', 'T cells', 'T cells',  # 6-10
                    'NA', 'NA', 'NA', 'NA', 'NA', 'NA'),   # 11-16
  'Resistant Fit' = c('T cells', 'Monocytes', 'T cells', 'B cells', 'T cells', 'Monocytes',   # 0-5
                      'NK cells', 'T cells', 'T cells', 'T cells', 'T cells',   # 6-10
                      'NA', 'Monocytes', 'B cells', 'NA', 'NA'),   # 11-15
  'SCRAN' = c('T cells', 'Monocytes', 'B cells', 'T cells', 'T cells', 'Monocytes',   # 0-5
              'NK cells', 'T cells', 'T cells', 'T cells', 'T cells',  # 6-10
              'NA', 'Monocytes', 'B cells', 'NA', 'NA', 'NA')   # 11-16
)

# Resistant Fit
# SCRAN

# -----------------------------
# Visulazation
# -----------------------------
# List of Title for each subgraph
title.list <- c('Log-scale', 'SCTransform', 'RF Norm', 'SCRAN')
g <- list()
label.list <- c('(a)', '(b)', '(c)', '(d)', '(e)')
for (i in seq(4)) {
  file_name <- paste0(OUT_DIR, file_list[i])
  res <- readRDS(file = file_name)
  ident <- celltype[[i]]
  cols <- c("T cells" = "#1B9E77", "B cells" = "#D95F02",
            "Monocytes" = "#7570B3", "NK cells" = "#E7298A")
  levels(res$md$seurat_clusters) <- ident
  res$md$cell <- rownames(res$md)
  select.barcode <- intersect(rownames(res$md), rownames(res$umap@cell.embeddings))
  p.data <- cbind(res$md[select.barcode, c('cell', 'seurat_clusters')], res$umap@cell.embeddings)
  p.data <- filter(p.data, seurat_clusters != 'NA')

  p <- ggplot()
  p <- p + geom_point(data = p.data,
                      aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters),
                      shape = 16, size = 0.2)
  p <- p + scale_colour_manual(values = cols,
                               aesthetics = c("colour", "fill"))
  p <- p + xlab('UMAP 1')
  p <- p + ylab('UMAP 2')
  p <- p + ggtitle(title.list[i])
  p <- p + theme(legend.position = "none")
  p <- p + theme(plot.title = element_text(size = 18,
                                           face = "bold",
                                           hjust = 0.5))
  #p <- p + theme(legend.title = element_blank())
  p <- p + labs(tag = label.list[i])
  p <- p + guides(colour = guide_legend(override.aes = list(size = 5)))
  p <- p + theme(axis.ticks=element_blank(), axis.text = element_blank())
  g[[length(g) + 1]] <- p
}
file_name <- paste0(OUT_DIR, 'scvi-embedding.csv')
p.data <- read.csv(file_name)
p.data <- p.data[!is.na(p.data$cluster),]
levels(p.data$cluster) <- c("B cells", "Monocytes", "NK cells", "T cells")
p <- ggplot()
p <- p + geom_point(data = p.data,
                    aes(x = umap_x, y = umap_y, color = cluster),
                    shape = 16, size = 0.15)
p <- p + scale_colour_manual(values = cols,
                             aesthetics = c("colour", "fill"))
p <- p + xlab('UMAP 1')
p <- p + ylab('UMAP 2')
p <- p + ggtitle('scVI')
p <- p + theme(legend.position = c(0.85, 0.2))
#p <- p + theme(legend.position = 'bottom')
p <- p + theme(plot.title = element_text(size = 18,
                                         face = "bold",
                                         hjust = 0.5))
p <- p + theme(legend.title = element_blank())
p <- p + labs(tag = label.list[5])
p <- p + guides(colour = guide_legend(override.aes = list(size = 5)))
p <- p + theme(axis.ticks=element_blank(), axis.text = element_blank()) 
#p
g[[length(g) + 1]] <- p

multiplot(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], cols = 5)