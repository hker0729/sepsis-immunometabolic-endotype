# scRNA increment: KLRF1+ vs KLRF1- NK cells in GSE167363
suppressMessages({
  library(Seurat)
})
set.seed(42)
OUT <- "H:/data/.openclaw/workspace/sepsis_project/output/03_immunosenescence/20_stitch/scRNA"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cat("Loading integrated object (330MB)...\n")
obj <- readRDS("H:/data/.openclaw/workspace/sepsis_project/output/07_single_cell/seurat_integrated.rds")
cat("Cells:", ncol(obj), "Genes:", nrow(obj), "\n")

# annotation
ann <- read.csv("H:/data/.openclaw/workspace/sepsis_project/output/07_single_cell/cluster_annotations.csv")
map <- setNames(ann$cell_type, as.character(ann$cluster))
obj$cell_type <- unname(map[as.character(obj$seurat_clusters)])
cat("cell_type table:\n"); print(table(obj$cell_type))

# NK subset
nk <- subset(obj, cell_type == "NK")
cat("NK cells:", ncol(nk), "\n")
cat("NK by sepsis:\n"); print(table(nk$sepsis))
cat("Joining layers...\n")
nk <- JoinLayers(nk, assay = "RNA")
cat("Layers after join:", paste(Layers(nk, assay = "RNA"), collapse = ", "), "\n")

# KLRF1 expression
get_expr <- function(obj, genes) {
  LayerData(obj, assay = "RNA", layer = "data")[genes, , drop = FALSE]
}
kl <- get_expr(nk, "KLRF1")[1, ]
cat("KLRF1 detected in NK:", sum(kl > 0), "/", length(kl), sprintf("(%.1f%%)\n", 100*mean(kl > 0)))
nk$klrf1_pos <- kl > 0

# gene sets
cytotox <- intersect(c("GZMB","PRF1","GNLY","NKG7","KLRD1","KLRK1","SPON2","CTSW","GZMA","GZMH"), rownames(nk))
oxphos <- grep("^(NDUF|SDH|UQCR|COX|ATP5)", rownames(nk), value = TRUE)
glyco <- intersect(c("HK1","HK2","HK3","PFKL","PFKM","PFKP","PKM","LDHA","LDHB","ENO1","ENO2","GAPDH","SLC2A1","SLC2A3","PGK1","PGM1","TPI1","ALDOA"), rownames(nk))
exh <- intersect(c("PDCD1","CTLA4","LAG3","TIGIT","HAVCR2"), rownames(nk))
tfs <- intersect(c("TBX21","EOMES"), rownames(nk))
cat(sprintf("gene sets: cytotox=%d oxphos=%d glyco=%d exh=%d tfs=%d\n",
            length(cytotox), length(oxphos), length(glyco), length(exh), length(tfs)))

expr_all <- get_expr(nk, unique(c(cytotox, oxphos, glyco, exh, tfs, "KLRF1")))
score_expr <- function(genes) {
  if (length(genes) < 3) return(rep(NA, ncol(nk)))
  colMeans(expr_all[genes, , drop = FALSE])
}
nk$cytotox_score <- score_expr(cytotox)
nk$oxphos_score <- score_expr(oxphos)
nk$glyco_score <- score_expr(glyco)
nk$exhaustion_score <- score_expr(exh)

md <- nk@meta.data
cat("\n=== KLRF1+ vs KLRF1- NK (all) ===\n")
for (sc in c("cytotox_score","oxphos_score","glyco_score","exhaustion_score")) {
  w <- wilcox.test(md[[sc]][md$klrf1_pos], md[[sc]][!md$klrf1_pos])
  cat(sprintf("%s: KLRF1+ mean=%.3f vs KLRF1- mean=%.3f, p=%.3e\n",
              sc, mean(md[[sc]][md$klrf1_pos]), mean(md[[sc]][!md$klrf1_pos]), w$p.value))
}
for (tf in tfs) {
  e <- expr_all[tf, ]
  w <- wilcox.test(e[md$klrf1_pos], e[!md$klrf1_pos])
  cat(sprintf("%s: KLRF1+ mean=%.3f vs KLRF1- mean=%.3f, p=%.3e\n",
              tf, mean(e[md$klrf1_pos]), mean(e[!md$klrf1_pos]), w$p.value))
}

cat("\n=== KLRF1+ NK proportion: Sepsis vs Control ===\n")
tab <- table(md$sepsis, md$klrf1_pos)
print(tab)
f <- fisher.test(tab)
cat(sprintf("Fisher p = %.3e; KLRF1+ proportion sepsis=%.1f%% vs control=%.1f%%\n",
            f$p.value, 100*tab["Sepsis","TRUE"]/sum(tab["Sepsis",]), 100*tab["Control","TRUE"]/sum(tab["Control",])))

cat("\n=== NK functional scores: Sepsis vs Control ===\n")
for (sc in c("cytotox_score","oxphos_score","glyco_score")) {
  w <- wilcox.test(md[[sc]][md$sepsis=="Sepsis"], md[[sc]][md$sepsis=="Control"])
  cat(sprintf("%s: Sepsis mean=%.3f vs Control mean=%.3f, p=%.3e\n",
              sc, mean(md[[sc]][md$sepsis=="Sepsis"]), mean(md[[sc]][md$sepsis=="Control"]), w$p.value))
}

# save
saveRDS(list(
  n_nk = ncol(nk),
  klrf1_pos_rate = mean(kl > 0),
  klrf1_pos_sepsis = as.numeric(tab["Sepsis","TRUE"]/sum(tab["Sepsis",])),
  klrf1_pos_control = as.numeric(tab["Control","TRUE"]/sum(tab["Control",])),
  fisher_p = f$p.value
), file.path(OUT, "klrf1_nk_summary.rds"))

# figures
png(file.path(OUT, "fig_klrf1_nk_scores.png"), width = 2600, height = 1800, res = 300)
par(mfrow = c(2, 2), mar = c(5, 5, 3, 1))
for (sc in c("cytotox_score","oxphos_score","glyco_score","exhaustion_score")) {
  boxplot(md[[sc]] ~ md$klrf1_pos, col = c("#95A5A6", "#C0392B"),
          names = c("KLRF1-", "KLRF1+"), outline = FALSE,
          ylab = sc, main = sprintf("%s by KLRF1 status (NK cells)", sc))
}
dev.off()

png(file.path(OUT, "fig_klrf1_nk_sepsis.png"), width = 2200, height = 1600, res = 300)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))
barplot(c(100*tab["Control","TRUE"]/sum(tab["Control",]), 100*tab["Sepsis","TRUE"]/sum(tab["Sepsis",])),
        names.arg = c("Control", "Sepsis"), col = c("#27AE60", "#C0392B"),
        ylab = "% KLRF1+ NK cells", main = sprintf("KLRF1+ NK proportion (Fisher p=%.1e)", f$p.value), ylim = c(0, max(100*tab/rowSums(tab))*1.3))
boxplot(md$cytotox_score ~ md$sepsis, col = c("#27AE60", "#C0392B"), names = c("Control","Sepsis"),
        outline = FALSE, ylab = "Cytotoxicity score", main = "NK cytotoxicity by sepsis")
dev.off()
cat("\nDone. Saved to", OUT, "\n")
