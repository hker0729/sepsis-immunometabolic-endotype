# C2. GSE95233: exprs already gene-symbol rownames
suppressMessages({library(pROC)})
OUT <- "H:/data/.openclaw/workspace/sepsis_project/ltie_manuscript/mps_preanalysis"
g <- readRDS("H:/data/.openclaw/workspace/sepsis_project/data/geo/GSE95233_processed_corrected.rds")
ex <- g$exprs; pd <- g$pdata
mps17 <- c("MTOR","BAX","BAK1","BID","RHOA","NLRP3","TNFRSF1A","TNFRSF1B",
           "RICTOR","RPTOR","MYD88","NFE2L2","TLR2","TLR4","TLR7","TLR8","NINJ1")
rn <- rownames(ex)
avail <- intersect(mps17, rn)
cat("genes available:", length(avail), "/17 ->", paste(avail, collapse=","), "\n")
cat("missing:", paste(setdiff(mps17, rn), collapse=","), "\n")
gmat <- ex[avail, , drop=FALSE]
# dedupe: if duplicate gene rows, take mean
if (any(duplicated(avail))) {
  cat("duplicated rows present, collapsing\n")
  gmat <- rowsum(gmat, group=avail, reorder=TRUE) / as.numeric(table(avail)[rownames(gmat)])
}
zscore_mean <- function(gm) { z <- t(scale(t(gm))); colMeans(z, na.rm=TRUE) }
pd$mps <- zscore_mean(gmat)
pd$dead <- ifelse(pd$survival == "Non Survivor", 1, ifelse(pd$survival == "Survivor", 0, NA))
d <- pd[!is.na(pd$dead), ]
cat("\nanalysis n =", nrow(d), "| deaths =", sum(d$dead), "| survivors =", sum(d$dead==0), "\n")
cat("group table:\n"); print(table(d$group))
m0 <- glm(dead ~ scale(mps), data=d, family=binomial)
cat(sprintf("MPS unadj OR=%.3f (%.3f-%.3f) p=%.4g\n", exp(coef(m0)[2]),
    exp(confint(m0)[2,1]), exp(confint(m0)[2,2]), summary(m0)$coefficients[2,4]))
m <- glm(dead ~ scale(mps) + age + gender, data=d, family=binomial)
cat(sprintf("MPS adj OR=%.3f (%.3f-%.3f) p=%.4g\n", exp(coef(m)[2]),
    exp(confint(m)[2,1]), exp(confint(m)[2,2]), summary(m)$coefficients[2,4]))
roc0 <- roc(d$dead, scale(d$mps), quiet=TRUE)
cat(sprintf("AUC=%.3f\n", as.numeric(auc(roc0))))
cat("\nsingle genes (adj age+sex):\n")
gmat_d <- gmat[, d$sample_id, drop=FALSE]
for (gene in c("NINJ1","BAK1","BAX","MYD88","TLR4","RICTOR","MTOR")) {
  if (!gene %in% rownames(gmat_d)) { cat(sprintf("  %-8s NA\n", gene)); next }
  m2 <- glm(dead ~ scale(gmat_d[gene, ]) + age + gender, data=d, family=binomial)
  cat(sprintf("  %-8s OR=%.3f p=%.4g\n", gene, exp(coef(m2)[2]), summary(m2)$coefficients[2,4]))
}
saveRDS(list(dat=d, gmat=gmat), file.path(OUT, "mps_GSE95233_scores.rds"))
cat("\nSaved.\n")
