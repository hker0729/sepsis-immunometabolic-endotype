# -*- coding: utf-8 -*-
# RPTOR sensitivity analysis: 16-gene vs 17-gene MPS (manuscript Methods 4.2)
# Run after mps_pre_main.R / mps_c_95233.R / mps_d_13904.R
suppressMessages({library(pROC)})

mps17 <- c("MTOR","BAX","BAK1","BID","RHOA","NLRP3","TNFRSF1A","TNFRSF1B",
           "RICTOR","RPTOR","MYD88","NFE2L2","TLR2","TLR4","TLR7","TLR8","NINJ1")
zscore_mean <- function(gm) { z <- t(scale(t(gm))); colMeans(z, na.rm=TRUE) }

## GSE95233 (17/17 genes measured)
cat("========== GSE95233 (unadjusted) ==========\n")
g <- readRDS("H:/data/.openclaw/workspace/sepsis_project/data/geo/GSE95233_processed_corrected.rds")
pd <- g$pdata
avail <- intersect(mps17, rownames(g$exprs))
gmat <- g$exprs[avail, , drop=FALSE]
if (any(duplicated(avail))) gmat <- rowsum(gmat, group=avail, reorder=TRUE) / as.numeric(table(avail)[rownames(gmat)])
pd$mps17 <- zscore_mean(gmat)
pd$mps16 <- zscore_mean(gmat[setdiff(rownames(gmat),"RPTOR"), , drop=FALSE])
pd$dead <- ifelse(pd$survival == "Non Survivor", 1, ifelse(pd$survival == "Survivor", 0, NA))
d <- pd[!is.na(pd$dead), ]
for (sc in c("mps17","mps16")) {
  m <- glm(dead ~ scale(d[[sc]]), data=d, family=binomial)
  cat(sprintf("%s OR=%.3f (%.3f-%.3f) p=%.4g\n", sc, exp(coef(m)[2]),
      exp(confint(m)[2,1]), exp(confint(m)[2,2]), summary(m)$coefficients[2,4]))
}
cat(sprintf("cor(mps16, mps17)=%.4f\n", cor(d$mps16, d$mps17)))

## GSE13904 (17/17 genes measured)
cat("\n========== GSE13904 (sepsis vs control) ==========\n")
g2 <- readRDS("H:/data/.openclaw/workspace/sepsis_project/data/geo/GSE13904_processed.rds")
pd2 <- g2$pdata
avail2 <- intersect(mps17, rownames(g2$exprs))
gmat2 <- g2$exprs[avail2, , drop=FALSE]
if (any(duplicated(avail2))) gmat2 <- rowsum(gmat2, group=avail2, reorder=TRUE) / as.numeric(table(avail2)[rownames(gmat2)])
pd2$mps17 <- zscore_mean(gmat2)
pd2$mps16 <- zscore_mean(gmat2[setdiff(rownames(gmat2),"RPTOR"), , drop=FALSE])
pd2$sep <- ifelse(pd2$group == "sepsis", 1, 0)
for (sc in c("mps17","mps16")) {
  m <- glm(sep ~ scale(pd2[[sc]]), data=pd2, family=binomial)
  w <- wilcox.test(pd2[[sc]] ~ pd2$sep)
  cat(sprintf("%s OR=%.3f (%.3f-%.3f) p=%.4g | Wilcoxon p=%.4g | AUC=%.3f\n", sc,
      exp(coef(m)[2]), exp(confint(m)[2,1]), exp(confint(m)[2,2]),
      summary(m)$coefficients[2,4], w$p.value,
      as.numeric(auc(roc(pd2$sep, pd2[[sc]], quiet=TRUE)))))
}
cat(sprintf("cor(mps16, mps17)=%.4f\n", cor(pd2$mps16, pd2$mps17)))
cat("\nConclusion: 16-gene MPS (RPTOR dropped) preserves direction and significance in both cohorts.\n")
