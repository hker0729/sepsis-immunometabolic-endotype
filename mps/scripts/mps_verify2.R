# Compute all remaining manuscript numbers from expression matrices
suppressMessages({library(survival)})
OUT <- "H:/data/.openclaw/workspace/sepsis_project/ltie_manuscript/mps_preanalysis"
mps17 <- c("MTOR","BAX","BAK1","BID","RHOA","NLRP3","TNFRSF1A","TNFRSF1B",
           "RICTOR","RPTOR","MYD88","NFE2L2","TLR2","TLR4","TLR7","TLR8","NINJ1")
imm_genes <- c("TLR2","TLR4","TLR7","TLR8","MYD88","NLRP3","TNFRSF1A","TNFRSF1B","NFE2L2")
dex_genes <- c("BAX","BAK1","BID","NINJ1","RHOA","RICTOR")
zscore_mean <- function(gm) { z <- t(scale(t(gm))); colMeans(z, na.rm=TRUE) }

cat("===== GSE65682: single-gene & module Cox =====\n")
g <- readRDS("H:/data/.openclaw/workspace/sepsis_project/data/geo/GSE65682_processed_full.rds")
ex <- g$exprs; pd <- g$pdata
avail <- intersect(mps17, rownames(ex))
cat("genes:", length(avail), "/17; missing:", paste(setdiff(mps17, rownames(ex)), collapse=","), "\n")
gmat <- ex[avail, , drop=FALSE]
if (any(duplicated(avail))) gmat <- rowsum(gmat, group=avail, reorder=TRUE) / as.numeric(table(avail)[rownames(gmat)])
# align samples
common <- intersect(colnames(gmat), pd$sample_id)
pd <- pd[match(common, pd$sample_id), ]
gmat <- gmat[, common]
pd$status <- as.numeric(pd$mortality_event_28days)
pd$time <- as.numeric(pd$time_to_event_28days)
cat("n =", nrow(pd), "deaths =", sum(pd$status), "\n")
for (gene in c("MYD88","TLR4","NINJ1","BAK1","BAX","RICTOR")) {
  m <- coxph(Surv(time, status) ~ scale(gmat[gene, ]), data=pd)
  s <- summary(m)
  cat(sprintf("  %-8s HR=%.3f (%.3f-%.3f) p=%.4g\n", gene, exp(coef(m)),
      exp(confint(m)[,1]), exp(confint(m)[,2]), s$coefficients[,"Pr(>|z|)"]))
}
pd$ImmuneSensing <- zscore_mean(gmat[intersect(imm_genes, rownames(gmat)), ])
pd$DeathExecution <- zscore_mean(gmat[intersect(dex_genes, rownames(gmat)), ])
r_mod <- cor(pd$ImmuneSensing, pd$DeathExecution, method="spearman")
cat(sprintf("  ImmuneSensing vs DeathExecution Spearman r = %.3f\n", r_mod))
for (mod in c("ImmuneSensing","DeathExecution")) {
  m <- coxph(Surv(time, status) ~ scale(pd[[mod]]), data=pd)
  s <- summary(m)
  cat(sprintf("  %-15s HR=%.3f (%.3f-%.3f) p=%.4g\n", mod, exp(coef(m)),
      exp(confint(m)[,1]), exp(confint(m)[,2]), s$coefficients[,"Pr(>|z|)"]))
}
# MPS from gmat for merge with immunity
pd$mps <- zscore_mean(gmat)
saveRDS(pd, file.path(OUT, "mps_GSE65682_full.rds"))

cat("\n===== MPS vs immunity correlation + IDI/NRI =====\n")
imm <- readRDS("H:/data/.openclaw/workspace/sepsis_project/output/03_immunosenescence/02_scoring/scores_GSE65682.rds")
d <- merge(pd[, c("sample_id","status","time","mps","ImmuneSensing","DeathExecution")],
           imm[, c("sample_id","immunity_score","age","sex")], by="sample_id")
r_mp <- cor(d$mps, d$immunity_score, method="spearman")
cat(sprintf("  MPS vs immunity_score Spearman r = %.3f\n", r_mp))
# KM tertiles
library(survival)
d$t3 <- cut(d$mps, breaks=quantile(d$mps, c(0,1/3,2/3,1)), include.lowest=TRUE, labels=c("T1","T2","T3"))
km <- survdiff(Surv(time, status) ~ t3, data=d)
cat(sprintf("  KM tertile log-rank p = %.4g\n", 1 - pchisq(km$chisq, df=2)))
# tertile death rates
print(round(tapply(d$status, d$t3, mean), 3))
# IDI/NRI via riskRegression or manual: use survIDINRI if available
has_survIDINRI <- requireNamespace("survIDINRI", quietly=TRUE)
cat("  survIDINRI available:", has_survIDINRI, "\n")
if (has_survIDINRI) {
  library(survIDINRI)
  t0 <- 28
  d2 <- d[complete.cases(d[, c("mps","immunity_score","age","sex","time","status")]), ]
  # base model: age+sex+immunity ; new: +MPS
  Xbase <- model.matrix(~ age + sex + immunity_score, data=d2)[, -1]
  Xnew  <- model.matrix(~ age + sex + immunity_score + mps, data=d2)[, -1]
  set.seed(42)
  res <- IDI.INF(indata=as.matrix(data.frame(time=d2$time, status=d2$status)),
                 covs0=Xbase, covs1=Xnew, t0=t0, npert=200)
  cat(sprintf("  IDI=%.4f (%.4f-%.4f) p=%.4g\n", res$IDI$Estimate, res$IDI$lower, res$IDI$upper, res$IDI$p))
  cat(sprintf("  NRI=%.4f (%.4f-%.4f) p=%.4g\n", res$NRI$Estimate, res$NRI$lower, res$NRI$upper, res$NRI$p))
}
