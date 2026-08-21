# D. GSE13904: MPS sepsis vs control
suppressMessages({library(pROC)})
OUT <- "H:/data/.openclaw/workspace/sepsis_project/ltie_manuscript/mps_preanalysis"
g <- readRDS("H:/data/.openclaw/workspace/sepsis_project/data/geo/GSE13904_processed.rds")
ex <- g$exprs; pd <- g$pdata
mps17 <- c("MTOR","BAX","BAK1","BID","RHOA","NLRP3","TNFRSF1A","TNFRSF1B",
           "RICTOR","RPTOR","MYD88","NFE2L2","TLR2","TLR4","TLR7","TLR8","NINJ1")
avail <- intersect(mps17, rownames(ex))
cat("genes available:", length(avail), "/17 ->", paste(avail, collapse=","), "\n")
cat("missing:", paste(setdiff(mps17, rownames(ex)), collapse=","), "\n")
gmat <- ex[avail, , drop=FALSE]
if (any(duplicated(avail))) {
  gmat <- rowsum(gmat, group=avail, reorder=TRUE) / as.numeric(table(avail)[rownames(gmat)])
}
zscore_mean <- function(gm) { z <- t(scale(t(gm))); colMeans(z, na.rm=TRUE) }
pd$mps <- zscore_mean(gmat)
pd$sep <- ifelse(pd$group == "sepsis", 1, 0)
w <- wilcox.test(pd$mps ~ pd$sep)
cat(sprintf("MPS sepsis (n=%d, mean=%.3f) vs control (n=%d, mean=%.3f): p=%.4g\n",
    sum(pd$sep==1), mean(pd$mps[pd$sep==1]), sum(pd$sep==0), mean(pd$mps[pd$sep==0]), w$p.value))
r <- roc(pd$sep, pd$mps, quiet=TRUE)
cat(sprintf("AUC=%.3f\n", as.numeric(auc(r))))
cat("\nsingle genes sepsis vs control (Wilcoxon):\n")
for (gene in c("NINJ1","BAK1","BAX","MYD88","TLR4","RICTOR","MTOR")) {
  if (!gene %in% rownames(gmat)) next
  v <- gmat[gene, ]
  w2 <- wilcox.test(v[pd$sep==1], v[pd$sep==0])
  cat(sprintf("  %-8s sepsis %.3f vs control %.3f | p=%.4g\n", gene,
      mean(v[pd$sep==1]), mean(v[pd$sep==0]), w2$p.value))
}
saveRDS(pd, file.path(OUT, "mps_GSE13904_scores.rds"))
cat("\nSaved.\n")
