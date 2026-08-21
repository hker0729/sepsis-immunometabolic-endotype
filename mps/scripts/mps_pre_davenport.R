# MPS validation in Davenport E-MTAB-4421 — z-score mean method (robust, NA-safe)
# consistent with GSE65682 pre-analysis rerun
suppressMessages(library(survival))
BASE <- "H:/data/.openclaw/workspace/sepsis_project/data/geo/E-MTAB-4421"
MAP  <- "H:/data/.openclaw/workspace/sepsis_project/ltie_manuscript/mps_preanalysis/mps17_probe_map_GPL6947.csv"

zscore_mean <- function(gmat) {
  z <- t(scale(t(gmat)))
  colMeans(z, na.rm=TRUE)
}

ex <- read.delim(file.path(BASE, "Davenport_sepsis_Jan2016_normalised_265.txt"), header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
sdrf <- read.delim(file.path(BASE, "E-MTAB-4421.sdrf.txt"), header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, sep="\t")
sdrf$id <- sdrf$`Comment[beadchiparray_id]`
sdrf$surv28 <- sdrf$`Characteristics[28 day survival]`
sdrf$age <- as.numeric(sdrf$`Characteristics[age]`)
sdrf$sex <- sdrf$`Characteristics[sex]`
common <- intersect(colnames(ex)[-1], sdrf$id)
sdrf <- sdrf[match(common, sdrf$id), ]
exm <- as.matrix(ex[, common]); rownames(exm) <- ex$Probe_ID

pmap <- read.csv(MAP, stringsAsFactors=FALSE)
gene_expr <- list()
for (i in seq_len(nrow(pmap))) {
  g <- pmap$gene[i]
  pr <- unlist(strsplit(pmap$probes[i], ","))
  pr <- pr[pr != "" & pr %in% rownames(exm)]
  if (length(pr) == 0) next
  gene_expr[[g]] <- colMeans(exm[pr, , drop=FALSE], na.rm=TRUE)
}
gmat <- do.call(rbind, gene_expr)
cat("genes:", nrow(gmat), "\n")

imm <- intersect(c("TLR2","TLR4","TLR7","TLR8","MYD88","NLRP3","TNFRSF1A","TNFRSF1B","NFE2L2"), rownames(gmat))
dex <- intersect(c("BAX","BAK1","BID","NINJ1","RHOA","RICTOR"), rownames(gmat))

dat <- data.frame(id=common,
                  surv28 = factor(ifelse(sdrf$surv28=="survivor", 0, 1)),
                  age=sdrf$age, sex=sdrf$sex,
                  MPS = zscore_mean(gmat),
                  ImmuneSensing = zscore_mean(gmat[imm, , drop=FALSE]),
                  DeathExecution = zscore_mean(gmat[dex, , drop=FALSE]))
dat <- dat[complete.cases(dat[, c("MPS","surv28")]), ]
cat("analysis n =", nrow(dat), "| deaths:", sum(dat$surv28==1), "| survivors:", sum(dat$surv28==0), "\n")

for (mod in c("MPS","ImmuneSensing","DeathExecution")) {
  m <- glm(surv28 ~ scale(dat[[mod]]), data=dat, family=binomial)
  ss <- summary(m); or <- exp(coef(m)[2]); ci <- exp(confint(m)[2,])
  cat(sprintf("%-15s OR=%.3f (%.3f-%.3f) p=%.3g\n", mod, or, ci[1], ci[2], ss$coefficients[2,4]))
}
m_adj <- glm(surv28 ~ scale(MPS) + age + sex, data=dat, family=binomial)
ss <- summary(m_adj); or <- exp(coef(m_adj)[2]); ci <- exp(confint(m_adj)[2,])
cat(sprintf("%-15s OR=%.3f (%.3f-%.3f) p=%.3g  [adj age+sex]\n", "MPS adj", or, ci[1], ci[2], ss$coefficients[2,4]))

cat("\nSingle-gene logistic:\n")
for (g in c("NINJ1","BAK1","MYD88","BAX","RICTOR","NLRP3","TLR4")) {
  if (!g %in% rownames(gmat)) { cat(sprintf("%-8s NA\n", g)); next }
  m <- glm(surv28 ~ scale(gmat[g, ]), data=dat, family=binomial)
  ss <- summary(m); or <- exp(coef(m)[2]); ci <- exp(confint(m)[2,])
  cat(sprintf("%-8s OR=%.3f (%.3f-%.3f) p=%.3g\n", g, or, ci[1], ci[2], ss$coefficients[2,4]))
}

saveRDS(dat, "H:/data/.openclaw/workspace/sepsis_project/ltie_manuscript/mps_preanalysis/mps_Davenport_scores.rds")
cat("\nSaved.\n")
