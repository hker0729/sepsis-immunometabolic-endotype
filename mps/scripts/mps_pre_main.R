# MPS pre-analysis: mitoxyperilysis 17-gene ssGSEA score vs 28-day mortality (GSE65682, n=479)
# Core 17 genes verified against Research Square preprint rs-9963383
suppressMessages({library(GSVA); library(survival)})

rds <- "H:/data/.openclaw/workspace/sepsis_project/data/geo/GSE65682_processed_full.rds"
x <- readRDS(rds)
exprs <- x$exprs; p <- x$pdata

# --- 17 core mitoxyperilysis genes (verified) ---
mps17 <- c("MTOR","BAX","BAK1","BID","RHOA","NLRP3","TNFRSF1A","TNFRSF1B",
           "RICTOR","RPTOR","MYD88","NFE2L2","TLR2","TLR4","TLR7","TLR8","NINJ1")

rn <- rownames(exprs)
cat("exprs rows:", length(rn), "| first:", paste(head(rn,3),collapse=","), "\n")
hit <- mps17[mps17 %in% rn]
miss <- setdiff(mps17, rn)
cat("gene match:", length(hit), "/17\n")
cat("matched:", paste(hit, collapse=","), "\n")
if (length(miss)>0) cat("MISSING:", paste(miss, collapse=","), "\n")

# survival data
dat <- data.frame(status = as.numeric(p$status), time = as.numeric(p$time), age = as.numeric(p$age))
cat("events:", sum(dat$status), "| n:", nrow(dat), "\n")

# --- ssGSEA score ---
gset <- list(Mitoxyperilysis = hit)
gsva_res <- tryCatch({
  gsva(exprs, gset, method="ssgsea", kcdf="Gaussian", verbose=FALSE)
}, error=function(e) {
  # new GSVA API
  param <- GSVA::ssgseaParam(exprData = exprs, geneSets = gset)
  gsva(param)
})
mps <- as.numeric(gsva_res[1, ])
dat$mps <- mps
cat("\nMPS summary:\n"); print(summary(mps))

# --- Association with 28-day mortality ---
# 1) continuous Cox per SD
sd_mps <- sd(mps)
cox1 <- coxph(Surv(time, status) ~ scale(mps), data=dat)
hr1 <- exp(coef(cox1)); ci1 <- exp(confint(cox1))
cat(sprintf("\n[1] Cox per SD: HR=%.3f (%.3f-%.3f) p=%.3g\n", hr1, ci1[1], ci1[2], summary(cox1)$coefficients[5]))

# 2) adjusted (age+sex)
dat$sex <- as.factor(p$gender)
cox2 <- coxph(Surv(time, status) ~ scale(mps) + age + sex, data=dat)
hr2 <- exp(coef(cox2)); ci2 <- exp(confint(cox2))
cat(sprintf("[2] Cox adj(age+sex): HR=%.3f (%.3f-%.3f) p=%.3g\n", hr2[1], ci2[1,1], ci2[1,2], summary(cox2)$coefficients[1,5]))

# 3) median split KM
med <- median(mps)
dat$grp <- ifelse(mps >= med, "HighMPS", "LowMPS")
km <- survdiff(Surv(time, status) ~ grp, data=dat)
p_km <- 1 - pchisq(km$chisq, 1)
cat(sprintf("[3] KM median split: events High=%d/%d Low=%d/%d, log-rank p=%.3g\n",
    sum(dat$status[dat$grp=="HighMPS"]), sum(dat$grp=="HighMPS"),
    sum(dat$status[dat$grp=="LowMPS"]), sum(dat$grp=="LowMPS"), p_km))

# 4) tertile trend
q3 <- quantile(mps, c(1/3, 2/3))
dat$t3 <- cut(mps, breaks=c(-Inf, q3[1], q3[2], Inf), labels=c("T1","T2","T3"))
cox_t <- coxph(Surv(time, status) ~ as.numeric(t3), data=dat)
cat(sprintf("[4] Tertile trend: HR per tertile=%.3f (%.3f-%.3f) p=%.3g\n",
    exp(coef(cox_t)), exp(confint(cox_t))[1], exp(confint(cox_t))[2], summary(cox_t)$coefficients[5]))
tb <- table(dat$t3, dat$status)
cat("tertile x status:\n"); print(tb)

# 5) correlation with key clinical (age)
cat(sprintf("\n[5] cor(MPS, age)=%.3f p=%.3g\n", cor.test(mps, dat$age)$estimate, cor.test(mps, dat$age)$p.value))

# --- save ---
outdir <- "H:/data/.openclaw/workspace/sepsis_project/ltie_manuscript/mps_preanalysis"
dir.create(outdir, showWarnings=FALSE)
saveRDS(data.frame(sample_id=p$sample_id, status=dat$status, time=dat$time, mps=mps, grp=dat$grp, t3=dat$t3),
        file.path(outdir, "mps_GSE65682_scores.rds"))
cat("\nSaved to", outdir, "\n")
