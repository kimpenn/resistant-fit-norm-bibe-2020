# Title     : Check overlap among markers
# Created by: derek
# Created on: 7/10/20

PRJ_DIR <- "/home/derek/research/Kim-Lab/normalization-simulation/"
IN_DIR <- paste0(PRJ_DIR, "data/pbmc-33k/filtered_gene_bc_matrices/hg19/")
OUT_DIR <- paste0(PRJ_DIR, "exp/exp-14/out/")

WRT_FLG <- TRUE

# -----------------------------
# Read data and Preprocess
# -----------------------------
file_name <- paste0(OUT_DIR, '2-marker.rf.cd4.cd8.eff.csv')
markers.rf <- read.csv(file = file_name, row.names=1)

file_name <- paste0(OUT_DIR, '2-marker.log.cd4.cd8.eff.csv')
markers.log <- read.csv(file = file_name, row.names=1)

file_name <- paste0(OUT_DIR, '2-marker.sct.cd4.cd8.eff.csv')
markers.sct <- read.csv(file = file_name, row.names=1)

file_name <- paste0(OUT_DIR, '2-marker.scran.cd4.cd8.eff.csv')
markers.scran <- read.csv(file = file_name, row.names=1)

file_name <- paste0(OUT_DIR, '2-marker.scvi.cd4.cd8.eff.csv')
markers.scvi <- read.csv(file = file_name, row.names=1)
# -----------------------------
# Read data and Preprocess
# -----------------------------
pos.markers.sct <- markers.sct[markers.sct$avg_logFC >0, ]
pos.markers.rf <- markers.rf[markers.rf$avg_logFC >0, ]
pos.markets.sct.rf <- intersect(pos.markers.sct$symbol, pos.markers.rf$symbol)

neg.markers.sct <- markers.sct[markers.sct$avg_logFC <0, ]
neg.markers.rf <- markers.rf[markers.rf$avg_logFC <0, ]
neg.markets.sct.rf <- intersect(neg.markers.sct$symbol, neg.markers.rf$symbol)
# -----------------------------
# Total Marker
# -----------------------------
markers.rf.log <- intersect(markers.rf$symbol, markers.log$symbol)
markers.rf.sct <- intersect(markers.rf$symbol, markers.sct$symbol)
markers.rf.scran <- intersect(markers.rf$symbol, markers.scran$symbol)
markers.rf.scvi <- intersect(markers.rf$symbol, markers.scvi$symbol)
length(markers.rf.log)
length(markers.rf.sct)
length(markers.rf.scran)
length(markers.rf.scvi)
# -----------------------------
# Neg Marker
# -----------------------------
neg.markers.rf.log <- intersect(markers.rf$symbol[markers.rf$avg_logFC < 0],
                            markers.log$symbol[markers.log$avg_logFC < 0])
neg.markers.rf.sct <- intersect(markers.rf$symbol[markers.rf$avg_logFC < 0],
                            markers.sct$symbol[markers.sct$avg_logFC < 0])
neg.markers.rf.scran <- intersect(markers.rf$symbol[markers.rf$avg_logFC < 0],
                              markers.scran$symbol[markers.scran$avg_logFC < 0])
neg.markers.rf.scvi <- intersect(markers.rf$symbol[markers.rf$avg_logFC < 0],
                              markers.scvi$symbol[markers.scvi$bayes_factor < 0])
length(neg.markers.rf.log)
length(neg.markers.rf.sct)
length(neg.markers.rf.scran)
length(neg.markers.rf.scvi)
setdiff(markers.scvi$symbol[markers.scvi$bayes_factor < 0],
        markers.rf$symbol[markers.rf$avg_logFC < 0])
# -----------------------------
# Pos Marker
# -----------------------------
pos.markers.rf.log <- intersect(markers.rf$symbol[markers.rf$avg_logFC > 0],
                            markers.log$symbol[markers.log$avg_logFC > 0])
pos.markers.rf.sct <- intersect(markers.rf$symbol[markers.rf$avg_logFC > 0],
                            markers.sct$symbol[markers.sct$avg_logFC > 0])
pos.markers.rf.scran <- intersect(markers.rf$symbol[markers.rf$avg_logFC > 0],
                              markers.scran$symbol[markers.scran$avg_logFC > 0])
pos.markers.rf.scvi <- intersect(markers.rf$symbol[markers.rf$avg_logFC > 0],
                              markers.scvi$symbol[markers.scvi$bayes_factor > 0])
length(pos.markers.rf.log)
length(pos.markers.rf.sct)
length(pos.markers.rf.scran)
length(pos.markers.rf.scvi)

