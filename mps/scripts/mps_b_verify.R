# Verify SRS interaction robustness
BASE <- "H:/data/.openclaw/workspace/sepsis_project/data/geo/E-MTAB-4421"
MAP  <- "H:/data/.openclaw/workspace/sepsis_project/ltie_manuscript/mps_preanalysis/mps17_probe_map_GPL6947.csv"
zscore_mean <- function(gmat) { z <- t(scale(t(gmat))); colMeans(z, na.rm=TRUE) }
ex <- read.delim(file.path(BASE, "Davenport_sepsis_Jan2016_normalised_265.txt"), header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
sdrf <- read.delim(file.path(BASE, "E-MTAB-4421.sdrf.txt"), header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, sep="\t")
sdrf$id <- sdrf[["Assay Name"]]
sdrf$srs <- sdrf[["Characteristics[sepsis response signature group]"]]
sdrf$surv28c <- sdrf[["Characteristics[28 day survival]"]]
sdrf$age <- as.numeric(sdrf[["Characteristics[age]"]])
sdrf$sex <- sdrf[["Characteristics[sex]"]]
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
d2 <- data.frame(id=common, srs2 = ifelse(sdrf$srs == "group 2", 1, 0),
                 dead = ifelse(sdrf$surv28c == "non-survivor", 1, 0),
                 age = sdrf$age, sex = sdrf$sex, MPS = zscore_mean(gmat))
d2 <- d2[complete.cases(d2[, c("MPS","dead","srs2")]), ]
cat("n =", nrow(d2), "| deaths =", sum(d2$dead), "\n")
cat("MPS by SRS: SRS1 mean", round(mean(d2$MPS[d2$srs2==0]),3), "| SRS2 mean", round(mean(d2$MPS[d2$srs2==1]),3), "\n")
w <- wilcox.test(d2$MPS[d2$srs2==0], d2$MPS[d2$srs2==1])
cat("MPS SRS1 vs SRS2 p =", format.pval(w$p.value), "\n")
cat("death rate: SRS1", round(mean(d2$dead[d2$srs2==0]),3), "| SRS2", round(mean(d2$dead[d2$srs2==1]),3), "\n")
cat("\n-- unadjusted models --\n")
m0 <- glm(dead ~ scale(MPS), data=d2, family=binomial)
cat(sprintf("MPS alone: OR=%.3f p=%.4g\n", exp(coef(m0)[2]), summary(m0)$coefficients[2,4]))
m1 <- glm(dead ~ scale(MPS) + srs2, data=d2, family=binomial)
cat(sprintf("MPS + SRS: OR=%.3f p=%.4g | SRS2 OR=%.3f p=%.4g\n",
    exp(coef(m1)[2]), summary(m1)$coefficients[2,4],
    exp(coef(m1)[3]), summary(m1)$coefficients[3,4]))
m2 <- glm(dead ~ scale(MPS) * srs2, data=d2, family=binomial)
ss <- summary(m2)$coefficients
cat(sprintf("MPS*SRS interaction (unadj): p=%.4g\n", ss[4,4]))
cat("  SRS1 MPS OR:", round(exp(ss[2,1]),3), "| SRS2 MPS OR:", round(exp(ss[2,1]+ss[4,1]),3), "\n")
# stratified unadjusted
for (srsv in 0:1) {
  dd <- d2[d2$srs2 == srsv, ]
  mm <- glm(dead ~ scale(MPS), data=dd, family=binomial)
  cat(sprintf("  SRS%d unadj: OR=%.3f p=%.4g\n", srsv+1, exp(coef(mm)[2]), summary(mm)$coefficients[2,4]))
}
cat("\n-- tertile of MPS by SRS --\n")
d2$t3 <- cut(d2$MPS, breaks=quantile(d2$MPS, c(0,1/3,2/3,1)), include.lowest=TRUE, labels=c("T1","T2","T3"))
print(table(d2$srs2, d2$t3))
print(round(tapply(d2$dead, list(d2$srs2, d2$t3), mean), 3))
