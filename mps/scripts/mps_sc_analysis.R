# MPS single-cell validation in GSE167363 (annotated object, 21k cells, 6 cell types)
suppressMessages({library(Seurat)})
s <- readRDS("H:/data/.openclaw/workspace/sepsis_project/output/03_immunosenescence/06_scRNA/seurat_scRNA.rds")

mps17 <- c("MTOR","BAX","BAK1","BID","RHOA","NLRP3","TNFRSF1A","TNFRSF1B",
           "RICTOR","RPTOR","MYD88","NFE2L2","TLR2","TLR4","TLR7","TLR8","NINJ1")
data <- LayerData(s, assay="RNA", layer="data")
g <- intersect(mps17, rownames(data))
cat("MPS genes available:", length(g), "/17 ->", paste(g, collapse=","), "\n")
cat("missing:", paste(setdiff(mps17, g), collapse=","), "\n\n")

# 1. Module score (full 17 where possible)
s <- AddModuleScore(s, features=list(g), name="MPS")
s$MPS <- s$MPS1

# 2. cell type assignment: argmax across module-score columns
ct_cols <- c("CD14_Mono1","NK_cell1","CD8_T_cell1","CD4_T_cell1","B_cell1","Neutrophil1")
score_mat <- as.matrix(s@meta.data[, ct_cols])
ct <- ct_cols[max.col(score_mat, ties.method="first")]
ct <- sub("1$", "", ct)
s$celltype <- factor(ct)
cat("celltype table:\n"); print(table(s$celltype))

# 3. Myeloid vs lymphoid MPS (Wilcoxon)
myeloid <- c("CD14_Mono","Neutrophil")
lymphoid <- c("NK_cell","CD8_T_cell","CD4_T_cell","B_cell")
my <- s$celltype %in% myeloid; ly <- s$celltype %in% lymphoid
w <- wilcox.test(s$MPS[my], s$MPS[ly])
cat(sprintf("\nMPS myeloid (n=%d, mean=%.3f) vs lymphoid (n=%d, mean=%.3f): W p=%.3g\n",
    sum(my), mean(s$MPS[my]), sum(ly), mean(s$MPS[ly]), w$p.value))

# 4. Per-celltype means
cat("\nPer-celltype MPS mean:\n")
print(round(tapply(s$MPS, s$celltype, mean), 3))

# 5. Sepsis vs Control within myeloid and overall
s$grp <- factor(ifelse(s$group=="Sepsis", "Sepsis", "Control"))
for (ctg in c("all","CD14_Mono","Neutrophil")) {
  idx <- if (ctg=="all") rep(TRUE, ncol(s)) else s$celltype == ctg
  if (sum(idx) < 10) next
  w <- wilcox.test(s$MPS[idx] ~ s$grp[idx])
  cat(sprintf("Sepsis vs Control MPS [%s]: n_sep=%d n_ctrl=%d | mean %.3f vs %.3f | p=%.3g\n",
      ctg, sum(s$grp[idx]=="Sepsis"), sum(s$grp[idx]=="Control"),
      mean(s$MPS[idx & s$grp=="Sepsis"]), mean(s$MPS[idx & s$grp=="Control"]), w$p.value))
}

# 6. Key single genes: myeloid vs lymphoid (log-normalized mean), and sepsis effect in myeloid
cat("\nSingle-gene mean expression by compartment (log-norm):\n")
for (gg in c("NINJ1","BAK1","BAX","MYD88","TLR4","MTOR","RICTOR")) {
  if (!gg %in% rownames(data)) { cat(sprintf("%-8s NA\n", gg)); next }
  v <- as.numeric(data[gg, ])
  my_v <- v[my]; ly_v <- v[ly]
  w <- wilcox.test(my_v, ly_v)
  cat(sprintf("%-8s myeloid %.3f vs lymphoid %.3f | p=%.3g\n", gg, mean(my_v), mean(ly_v), w$p.value))
}
cat("\nSingle-gene: sepsis vs control in MYELOID:\n")
for (gg in c("NINJ1","BAK1","BAX","MYD88","TLR4")) {
  if (!gg %in% rownames(data)) next
  v <- as.numeric(data[gg, ])
  w <- wilcox.test(v[my & s$grp=="Sepsis"], v[my & s$grp=="Control"])
  cat(sprintf("%-8s sepsis %.3f vs control %.3f | p=%.3g\n", gg,
      mean(v[my & s$grp=="Sepsis"]), mean(v[my & s$grp=="Control"]), w$p.value))
}

saveRDS(s, "H:/data/.openclaw/workspace/sepsis_project/ltie_manuscript/mps_preanalysis/mps_sc_GSE167363.rds")
cat("\nSaved.\n")
