# ============================================================================
# 04_immune_landscape_v2_quadrant.R - 免疫景观 (A 路线四象限)
# 2026-08-13 | 修复: 表达矩阵 = 真 GSE65682 (data/geo/GSE65682_processed_full.rds)
# 分组: Q1-Q4 四象限; 比较: Q1 (免疫代谢衰竭) vs 其他 / KW 四组
# ============================================================================
PROJECT <- "H:/data/.openclaw/workspace/sepsis_project"
setwd(PROJECT)

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(reshape2)

OUT <- "output/03_immunosenescence/04_immune"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat(" STAGE 4 v2: Immune Landscape (quadrant)\n")
cat("========================================\n\n")

scores <- readRDS("output/03_immunosenescence/02_scoring/scores_GSE65682.rds")
gene_sets <- readRDS("output/03_immunosenescence/01_gene_sets/immunosenescence_gene_sets.rds")

# 真 GSE65682 表达矩阵
dat <- readRDS("data/geo/GSE65682_processed_full.rds")
exprs <- dat$exprs
cat(sprintf("Expression: %d genes x %d samples\n", nrow(exprs), ncol(exprs)))

# 四象限列
if (!("quad_short" %in% colnames(scores))) {
  met_med <- median(scores$metabolism_score_pc1, na.rm = TRUE)
  imm_med <- median(scores$immunity_score, na.rm = TRUE)
  scores$met_hi <- ifelse(scores$metabolism_score_pc1 > met_med, 1, 0)
  scores$imm_hi <- ifelse(scores$immunity_score > imm_med, 1, 0)
  scores$quad_short <- factor(paste0(scores$met_hi, scores$imm_hi),
                              levels = c("00", "10", "01", "11"),
                              labels = c("Q1", "Q2", "Q3", "Q4"))
}

quad_cols <- c("Q1" = "#C0392B", "Q2" = "#E67E22", "Q3" = "#2E86C1", "Q4" = "#27AE60")

# ===========================================================================
# 1. 免疫细胞签名 (marker mean expression)
# ===========================================================================
cat("[1] Immune cell signatures...\n")
immune_signatures <- list(
  NK_cells = intersect(gene_sets$immunity$nk_effector, rownames(exprs)),
  CD8_T_cells = intersect(c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1"), rownames(exprs)),
  CD4_T_cells = intersect(c("CD4", "CD3E", "CD3D", "CD3G", "IL7R"), rownames(exprs)),
  B_cells = intersect(c("CD19", "MS4A1", "CD79A", "CD79B"), rownames(exprs)),
  Monocytes = intersect(c("CD14", "FCGR3A", "CSF1R", "CD163"), rownames(exprs)),
  Neutrophils = intersect(c("FCGR3B", "CSF3R", "CXCR2", "CEACAM8"), rownames(exprs)),
  Tregs = intersect(c("FOXP3", "IL2RA", "CTLA4", "IKZF2"), rownames(exprs)),
  Exhausted_T = intersect(c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4"), rownames(exprs))
)

ciber_res <- matrix(NA, ncol(exprs), length(immune_signatures))
colnames(ciber_res) <- names(immune_signatures)
rownames(ciber_res) <- colnames(exprs)
for (sig_name in names(immune_signatures)) {
  sig_genes <- immune_signatures[[sig_name]]
  if (length(sig_genes) >= 3) {
    ciber_res[, sig_name] <- colMeans(exprs[sig_genes, , drop = FALSE], na.rm = TRUE)
  }
}
cat(sprintf("  %d samples x %d cell types\n", nrow(ciber_res), ncol(ciber_res)))

# ===========================================================================
# 2. 免疫细胞差异 (KW 四象限 + Q1 vs others Wilcoxon)
# ===========================================================================
common <- intersect(rownames(ciber_res), rownames(scores))
cb_sub <- ciber_res[common, , drop = FALSE]
sc_sub <- scores[common, ]
sc_sub$quad_short <- factor(sc_sub$quad_short)

diff_res <- data.frame()
for (cell in colnames(cb_sub)) {
  kw <- kruskal.test(cb_sub[, cell] ~ sc_sub$quad_short)
  q1 <- sc_sub$quad_short == "Q1"
  wt <- wilcox.test(cb_sub[q1, cell], cb_sub[!q1, cell])
  meds <- tapply(cb_sub[, cell], sc_sub$quad_short, median, na.rm = TRUE)
  diff_res <- rbind(diff_res, data.frame(
    cell_type = cell,
    med_Q1 = meds["Q1"], med_Q2 = meds["Q2"], med_Q3 = meds["Q3"], med_Q4 = meds["Q4"],
    KW_p = kw$p.value, Q1_vs_others_p = wt$p.value
  ))
}
diff_res$KW_FDR <- p.adjust(diff_res$KW_p, method = "BH")
diff_res$Q1_FDR <- p.adjust(diff_res$Q1_vs_others_p, method = "BH")
cat("\n=== Immune cell differences by quadrant ===\n")
print(diff_res, digits = 3, row.names = FALSE)
write.csv(diff_res, file.path(OUT, "immune_cells_by_quadrant.csv"), row.names = FALSE)

# 热图 (z-score; 行=样本, 列=细胞类型)
hm_data <- t(scale(t(cb_sub)))
hm_ord <- order(sc_sub$quad_short)
annot <- data.frame(Quadrant = sc_sub$quad_short, row.names = rownames(sc_sub))
annot_colors <- list(Quadrant = quad_cols)
pdf(file.path(OUT, "immune_cells_heatmap_quadrant.pdf"), width = 8, height = 6)
pheatmap(hm_data[hm_ord, ], annotation_row = annot,
         annotation_colors = annot_colors, show_rownames = FALSE,
         cluster_rows = FALSE, main = "Immune cell signatures by quadrant (GSE65682)",
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100))
dev.off()
cat("Saved: immune_cells_heatmap_quadrant.pdf\n")

# ===========================================================================
# 3. SASP + checkpoint 按象限
# ===========================================================================
cat("\n=== SASP & checkpoint by quadrant ===\n")
sasp_genes <- intersect(gene_sets$senescence$sasp, rownames(exprs))
ckpt_genes <- intersect(c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "BTLA",
                          "CD274", "PDCD1LG2", "VSIR"), rownames(exprs))
sasp_mat <- colMeans(exprs[sasp_genes, common, drop = FALSE], na.rm = TRUE)
ckpt_mat <- colMeans(exprs[ckpt_genes, common, drop = FALSE], na.rm = TRUE)

for (nm in c("sasp_mat", "ckpt_mat")) {
  v <- get(nm)
  kw <- kruskal.test(v ~ sc_sub$quad_short)
  meds <- tapply(v, sc_sub$quad_short, median, na.rm = TRUE)
  cat(sprintf("%s: KW p=%.4f | medians Q1-Q4: %.3f, %.3f, %.3f, %.3f\n",
              nm, kw$p.value, meds["Q1"], meds["Q2"], meds["Q3"], meds["Q4"]))
}

# 箱线图: 免疫细胞 (Q1 vs Q3+Q4 高免组)
plot_df <- melt(as.data.frame(cb_sub[, c("NK_cells", "CD8_T_cells", "CD4_T_cells", "Exhausted_T")]),
                id.vars = NULL, variable.name = "cell_type", value.name = "value")
plot_df$quad <- sc_sub$quad_short[rownames(cb_sub)]
p <- ggplot(plot_df, aes(x = quad, y = value, fill = quad)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.8) +
  facet_wrap(~cell_type, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = quad_cols) +
  labs(x = NULL, y = "Signature score (mean expr)", title = "Key immune signatures by quadrant") +
  theme_minimal(base_size = 10) + theme(legend.position = "none",
                                        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT, "immune_cells_boxplot_quadrant.pdf"), p, width = 10, height = 4)
cat("Saved: immune_cells_boxplot_quadrant.pdf\n")

cat("\n=== STAGE 4 v2 complete ===\n")
