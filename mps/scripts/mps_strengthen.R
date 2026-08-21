# MPS strengthen analyses: A) GSE65682 immune-axis adjustment  B) Davenport SRS interaction
# C) GSE95233 third prognostic cohort  D) GSE13904 sepsis-vs-control MPS
suppressMessages({library(survival); library(pROC)})
OUT <- "H:/data/.openclaw/workspace/sepsis_project/ltie_manuscript/mps_preanalysis"

## ---------- A. GSE65682: MPS independent of immunity_score? ----------
cat("========== A. GSE65682 immune-axis adjustment ==========\n")
mps <- readRDS(file.path(OUT, "mps_GSE65682_scores.rds"))   # sample_id, status, time, mps, grp, t3
imm <- readRDS("H:/data/.openclaw/workspace/sepsis_project/output/03_immunosenescence/02_scoring/scores_GSE65682.rds")
d <- merge(mps, imm[, c("sample_id","immunity_score","quadrant","quad_short","age","sex")], by="sample_id")
d$mps_sd <- scale(d$mps)[,1]; d$imm_sd <- scale(d$immunity_score)[,1]
d$status <- as.numeric(d$status); d$time <- as.numeric(d$time)
cat("n =", nrow(d), "| deaths =", sum(d$status==1), "\n")
m1 <- coxph(Surv(time, status) ~ mps_sd + age + sex, data=d)
m2 <- coxph(Surv(time, status) ~ mps_sd + imm_sd + age + sex, data=d)
m3 <- coxph(Surv(time, status) ~ mps_sd + imm_sd + mps_sd:imm_sd + age + sex, data=d)
cat("M1 MPS alone (adj age sex): HR=%.3f (%.3f-%.3f) p=%.4g\n")
s <- summary(m1); cat(sprintf("  MPS HR=%.3f (%.3f-%.3f) p=%.4g\n", exp(coef(m1)["mps_sd"]), exp(confint(m1)["mps_sd",1]), exp(confint(m1)["mps_sd",2]), s$coefficients["mps_sd","Pr(>|z|)"]))
s <- summary(m2); cat(sprintf("M2 + immunity: MPS HR=%.3f (%.3f-%.3f) p=%.4g | IMM HR=%.3f p=%.4g\n",
    exp(coef(m2)["mps_sd"]), exp(confint(m2)["mps_sd",1]), exp(confint(m2)["mps_sd",2]), s$coefficients["mps_sd","Pr(>|z|)"],
    exp(coef(m2)["imm_sd"]), s$coefficients["imm_sd","Pr(>|z|)"]))
s <- summary(m3); cat(sprintf("M3 + interaction: MPSxIMM p=%.4g\n", s$coefficients["mps_sd:imm_sd","Pr(>|z|)"]))
# quadrant-stratified MPS effect
for (q in sort(unique(d$quad_short))) {
  dd <- d[d$quad_short == q, ]
  if (sum(dd$status==1) < 5 || nrow(dd) < 30) next
  mm <- coxph(Surv(time, status) ~ mps_sd, data=dd)
  ss <- summary(mm)
  cat(sprintf("  MPS in %s (n=%d, deaths=%d): HR=%.3f p=%.4g\n", q, nrow(dd), sum(dd$status==1), exp(coef(mm)), ss$coefficients[,"Pr(>|z|)"]))
}

## ---------- B. Davenport: SRS interaction ----------
cat("\n========== B. Davenport SRS1/2 interaction ==========\n")
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
imm <- intersect(c("TLR2","TLR4","TLR7","TLR8","MYD88","NLRP3","TNFRSF1A","TNFRSF1B","NFE2L2"), rownames(gmat))
dex <- intersect(c("BAX","BAK1","BID","NINJ1","RHOA","RICTOR"), rownames(gmat))
d2 <- data.frame(id=common, srs2 = ifelse(sdrf$srs == "group 2", 1, 0),
                 dead = ifelse(sdrf$surv28c == "non-survivor", 1, 0),
                 age = sdrf$age, sex = sdrf$sex,
                 MPS = zscore_mean(gmat),
                 ImmuneSensing = zscore_mean(gmat[imm, , drop=FALSE]),
                 DeathExecution = zscore_mean(gmat[dex, , drop=FALSE]))
for (g in c("BAX","BAK1","NINJ1","MYD88")) if (g %in% rownames(gmat)) d2[[g]] <- gmat[g, ]
d2 <- d2[complete.cases(d2[, c("MPS","dead","srs2")]), ]
cat("n =", nrow(d2), "| deaths =", sum(d2$dead), "| SRS1 =", sum(d2$srs2==0), "| SRS2 =", sum(d2$srs2==1), "\n")
for (v in c("MPS","BAX","BAK1","NINJ1","MYD88","ImmuneSensing","DeathExecution")) {
  d2[[paste0(v,"_sd")]] <- scale(d2[[v]])[,1]
}
for (v in c("MPS","BAX","BAK1","NINJ1","MYD88")) {
  f <- as.formula(sprintf("dead ~ %s_sd + srs2 + %s_sd:srs2 + age + sex", v, v))
  mm <- glm(f, data=d2, family=binomial)
  ss <- summary(mm)$coefficients
  cat(sprintf("%s main OR=%.3f p=%.4g | SRS interaction p=%.4g\n", v,
      exp(coef(mm)[2]), ss[2,4], ss[4,4]))
  for (srsv in 0:1) {
    dd <- d2[d2$srs2 == srsv, ]
    mm2 <- glm(as.formula(sprintf("dead ~ %s_sd + age + sex", v)), data=dd, family=binomial)
    ci <- tryCatch(exp(confint(mm2)[2,]), error=function(e) c(NA,NA))
    cat(sprintf("   in SRS%d: OR=%.3f (%.3f-%.3f) p=%.4g\n", srsv+1,
        exp(coef(mm2)[2]), ci[1], ci[2], summary(mm2)$coefficients[2,4]))
  }
}
saveRDS(d2, file.path(OUT, "mps_Davenport_srs.rds"))

## ---------- C. GSE95233 third cohort ----------
cat("\n========== C. GSE95233 (Agilent GPL4204) ==========\n")
g <- readRDS("H:/data/.openclaw/workspace/sepsis_project/data/geo/GSE95233_processed_corrected.rds")
pd <- g$pdata
cat("pdata dim:", dim(pd), "| is_sepsis:", sum(pd$is_sepsis), "| is_healthy:", sum(pd$is_healthy), "\n")
print(table(pd$corrected_survival, useNA="ifany"))
print(table(pd$group, useNA="ifany"))
cat("cols:", paste(colnames(pd), collapse=" | "), "\n")
# probe annotation
annot <- read.delim(gzfile("H:/data/.openclaw/workspace/sepsis_project/data/geo/GPL4204.annot.gz"), stringsAsFactors=FALSE, check.names=FALSE, comment.char="#")
cat("GPL4204 annot dim:", dim(annot), "| cols:", paste(colnames(annot), collapse=" | "), "\n")
