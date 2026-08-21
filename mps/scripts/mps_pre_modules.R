# MPS module decomposition: immune-sensing vs death-execution modules
suppressMessages({library(GSVA); library(survival)})
x <- readRDS("H:/data/.openclaw/workspace/sepsis_project/data/geo/GSE65682_processed_full.rds")
exprs <- x$exprs; p <- x$pdata
dat <- data.frame(status=as.numeric(p$status), time=as.numeric(p$time), age=as.numeric(p$age), sex=as.factor(p$gender))

imm_sensing <- c("TLR2","TLR4","TLR7","TLR8","MYD88","NLRP3","TNFRSF1A","TNFRSF1B","NFE2L2")  # 9
death_exec  <- c("BAX","BAK1","BID","NINJ1","RHOA","MTOR","RICTOR")                            # 7
nij1_only   <- "NINJ1"
bax_family  <- c("BAX","BAK1","BID")

gs <- list(ImmuneSensing=imm_sensing[imm_sensing %in% rownames(exprs)],
           DeathExecution=death_exec[death_exec %in% rownames(exprs)],
           BAXfamily=bax_family[bax_family %in% rownames(exprs)],
           NINJ1=nij1_only[nij1_only %in% rownames(exprs)])
cat("module sizes:", sapply(gs, length), "\n")

sc <- tryCatch(gsva(exprs, gs, method="ssgsea", kcdf="Gaussian", verbose=FALSE),
               error=function(e) gsva(GSVA::ssgseaParam(exprData=exprs, geneSets=gs)))
for (mod in names(gs)) {
  s <- as.numeric(sc[mod, ])
  cox <- coxph(Surv(time, status) ~ scale(s), data=dat)
  hr <- exp(coef(cox)); ci <- exp(confint(cox))
  cat(sprintf("%-15s HR=%.3f (%.3f-%.3f) p=%.3g\n", mod, hr, ci[1], ci[2], summary(cox)$coefficients[5]))
}
# single genes Cox (key ones)
cat("\nSingle-gene Cox:\n")
for (g in c("NINJ1","BAX","BAK1","BID","TLR4","MYD88","NLRP3","MTOR","RICTOR")) {
  if (!g %in% rownames(exprs)) { cat(sprintf("%-8s NA\n", g)); next }
  s <- as.numeric(exprs[g, ])
  cox <- coxph(Surv(time, status) ~ scale(s), data=dat)
  hr <- exp(coef(cox)); ci <- exp(confint(cox))
  cat(sprintf("%-8s HR=%.3f (%.3f-%.3f) p=%.3g\n", g, hr, ci[1], ci[2], summary(cox)$coefficients[5]))
}
# correlation between modules
cat("\ncor(ImmuneSensing, DeathExecution):", cor(as.numeric(sc["ImmuneSensing",]), as.numeric(sc["DeathExecution",])), "\n")
