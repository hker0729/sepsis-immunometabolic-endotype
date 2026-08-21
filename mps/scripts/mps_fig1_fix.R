# Fig1 full regenerate with VERIFIED E-panel numbers (fixes hardcoded 0.80/0.012 & Davenport CI)
suppressMessages({library(survival); library(pROC)})
OUT <- "H:/data/.openclaw/workspace/sepsis_project/ltie_manuscript/mps_preanalysis"
FIG <- "H:/data/.openclaw/workspace/sepsis_project/ltie_manuscript/mps_manuscript/figures"

g65682 <- readRDS(file.path(OUT, "mps_GSE65682_full.rds"))
imm <- readRDS("H:/data/.openclaw/workspace/sepsis_project/output/03_immunosenescence/02_scoring/scores_GSE65682.rds")
d65682 <- merge(g65682[, c("sample_id","status","time","mps","ImmuneSensing","DeathExecution")],
                imm[, c("sample_id","immunity_score","age","sex")], by="sample_id")
g95233 <- readRDS(file.path(OUT, "mps_GSE95233_scores.rds"))$dat

pdf(file.path(FIG, "Fig1_MPS_landscape.pdf"), width=11, height=8)
par(mfrow=c(2,3), mar=c(4.5,4.5,3,1))

## A: schematic of the 17-gene program
plot(NA, xlim=c(0,10), ylim=c(0,10), axes=FALSE, xlab="", ylab="", main="Mitoxyperilysis program (17 genes)")
rect(0.3, 7.6, 4.2, 9.6, col=rgb(0.3,0.5,0.9,0.15), border="steelblue", lwd=2)
text(2.25, 9.2, "Immune sensing (9)", font=2, cex=0.9, col="steelblue")
text(2.25, 8.6, "TLR2 TLR4 TLR7 TLR8\nMYD88 NLRP3\nTNFRSF1A TNFRSF1B NFE2L2", cex=0.75)
rect(0.3, 0.6, 4.2, 4.0, col=rgb(0.9,0.3,0.3,0.12), border="firebrick", lwd=2)
text(2.25, 3.6, "Death execution (8)", font=2, cex=0.9, col="firebrick")
text(2.25, 3.0, "BAX BAK1 BID NINJ1\nRHOA RICTOR MTOR RPTOR", cex=0.75)
arrows(2.25, 7.6, 2.25, 4.0, lwd=2, col="grey40")
text(2.6, 5.8, "TLR/MYD88 signal", cex=0.7, col="grey30")
rect(5.0, 5.2, 9.6, 6.6, col=rgb(0.6,0.6,0.6,0.1), border="grey50", lwd=1.5)
text(7.3, 6.25, "Mitochondrial oxidative stress", cex=0.75)
text(7.3, 5.55, "BAX/BAK1/BID -> MOMP", cex=0.7)
rect(5.0, 3.4, 9.6, 4.8, col=rgb(0.6,0.6,0.6,0.1), border="grey50", lwd=1.5)
text(7.3, 4.45, "Sustained membrane contact", cex=0.75)
text(7.3, 3.75, "local oxidative damage", cex=0.7)
rect(5.0, 1.6, 9.6, 3.0, col=rgb(0.9,0.3,0.3,0.12), border="firebrick", lwd=1.5)
text(7.3, 2.65, "NINJ1-mediated", cex=0.75)
text(7.3, 1.95, "membrane rupture", cex=0.75, font=2)
arrows(7.3, 6.6, 7.3, 4.8, lwd=1.5, col="grey50"); arrows(7.3, 4.8, 7.3, 3.0, lwd=1.5, col="grey50")
text(7.3, 8.2, "mTORC2 regulation (RICTOR/MTOR/RPTOR)", cex=0.7, col="grey30")
arrows(7.3, 8.0, 7.3, 6.6, lwd=1.2, lty=2, col="grey40")

## B: GSE65682 forest (dynamic)
m <- coxph(Surv(time, status) ~ scale(mps), data=d65682)
mi <- coxph(Surv(time, status) ~ scale(ImmuneSensing), data=d65682)
md <- coxph(Surv(time, status) ~ scale(DeathExecution), data=d65682)
genes <- c("MYD88","BAK1","NINJ1","BAX","TLR4")
gmat <- readRDS("H:/data/.openclaw/workspace/sepsis_project/data/geo/GSE65682_processed_full.rds")$exprs
ghr <- gci <- gp <- c()
for (gn in genes) {
  mm <- coxph(Surv(time, status) ~ scale(gmat[gn, d65682$sample_id]), data=d65682)
  ghr <- c(ghr, exp(coef(mm))); gci <- rbind(gci, exp(confint(mm))); gp <- c(gp, summary(mm)$coefficients[,"Pr(>|z|)"])
}
labels <- c("MPS","Immune sensing","Death execution","MYD88","BAK1","NINJ1","BAX","TLR4")
hrs <- c(exp(coef(m)), exp(coef(mi)), exp(coef(md)), ghr)
los <- c(exp(confint(m))[1], exp(confint(mi))[1], exp(confint(md))[1], gci[,1])
his <- c(exp(confint(m))[2], exp(confint(mi))[2], exp(confint(md))[2], gci[,2])
ps  <- c(summary(m)$coefficients["scale(mps)","Pr(>|z|)"], summary(mi)$coefficients[,"Pr(>|z|)"], summary(md)$coefficients[,"Pr(>|z|)"], gp)
ord <- rev(order(hrs))
plot(NA, xlim=c(0.4,1.6), ylim=c(0.5, length(labels)+0.5), xlab="HR per SD (28-day mortality)", ylab="", yaxt="n", main="GSE65682 (n=479, 114 deaths)")
abline(v=1, lty=2, col="grey50")
for (i in seq_along(labels)) {
  j <- ord[i]; yy <- length(labels) - i + 1
  segments(los[j], yy, his[j], yy, lwd=2, col=ifelse(hrs[j]<1,"steelblue","firebrick"))
  points(hrs[j], yy, pch=18, cex=1.4, col=ifelse(hrs[j]<1,"steelblue","firebrick"))
  text(par("usr")[1]+0.02, yy, sprintf("%s  %.2f (%.2f-%.2f) p=%.3g", labels[j], hrs[j], los[j], his[j], ps[j]), adj=0, cex=0.8)
}

## C: KM tertiles
d65682$t3 <- cut(d65682$mps, breaks=quantile(d65682$mps, c(0,1/3,2/3,1)), include.lowest=TRUE, labels=c("T1 (low)","T2","T3 (high)"))
fit <- survfit(Surv(time, status) ~ t3, data=d65682)
plot(fit, col=c("firebrick","grey50","steelblue"), lwd=2, xlab="Days", ylab="Survival", main="GSE65682 by MPS tertile")
legend("bottomleft", legend=c("T1 (low MPS)","T2","T3 (high MPS)"), col=c("firebrick","grey50","steelblue"), lwd=2, bty="n", cex=0.8)
text(2, 0.15, "log-rank p = 0.0013", adj=0)

## D: GSE95233 ROC
r <- roc(g95233$dead, g95233$mps, quiet=TRUE)
plot(r, main=sprintf("GSE95233 (n=102, 34 deaths)\nAUC=%.3f", as.numeric(auc(r))))
mtext("MPS unadjusted OR=0.449, p=0.001", side=3, line=-1.5, cex=0.8)

## E: cross-cohort (VERIFIED numbers)
cc <- data.frame(cohort=c("GSE65682","GSE95233","Davenport"),
                 est=c(0.772, 0.449, 1.380), lo=c(0.651, 0.268, 0.994), hi=c(0.915, 0.708, 1.981),
                 metric=c("HR per SD","OR per SD","OR per SD"), p=c(0.003, 0.001, 0.068))
plot(NA, xlim=c(0.2,2.2), ylim=c(0.5,3.5), xlab="HR/OR per SD (28-day mortality)", ylab="", yaxt="n", main="Cross-cohort MPS effect")
abline(v=1, lty=2, col="grey50")
for (i in 1:3) {
  yy <- 3.5 - i
  segments(cc$lo[i], yy, cc$hi[i], yy, lwd=3, col=ifelse(cc$est[i]<1,"steelblue","firebrick"))
  points(cc$est[i], yy, pch=18, cex=2, col=ifelse(cc$est[i]<1,"steelblue","firebrick"))
  text(cc$lo[i]-0.05, yy, sprintf("%s\n%s=%.3f (%.2f-%.2f), p=%.3g", cc$cohort[i], cc$metric[i], cc$est[i], cc$lo[i], cc$hi[i], cc$p[i]), adj=1, cex=0.8)
}
text(1.15, 0.8, "Davenport: divergent (endotype-dependent, see Fig 5)", cex=0.7, col="firebrick")

## F: module HRs
plot(NA, xlim=c(0.5,1.2), ylim=c(0.5,2.5), xlab="HR per SD", ylab="", yaxt="n", main="Module decomposition (GSE65682)")
abline(v=1, lty=2, col="grey50")
mds <- list(ImmuneSensing=mi, DeathExecution=md)
for (i in seq_along(mds)) {
  yy <- 2.5 - i
  mm <- mds[[i]]
  segments(exp(confint(mm))[1], yy, exp(confint(mm))[2], yy, lwd=3, col="steelblue")
  points(exp(coef(mm)), yy, pch=18, cex=2, col="steelblue")
  text(0.52, yy, sprintf("%s  HR=%.3f, p=%.4g", names(mds)[i], exp(coef(mm)), summary(mm)$coefficients[,"Pr(>|z|)"]), adj=0)
}
dev.off()
cat("Fig1 regenerated with verified E panel\n")
