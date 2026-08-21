# Regenerate Fig1E, Fig2, Fig3 with verified numbers
suppressMessages({library(survival); library(pROC)})
OUT <- "H:/data/.openclaw/workspace/sepsis_project/ltie_manuscript/mps_preanalysis"
FIG <- "H:/data/.openclaw/workspace/sepsis_project/ltie_manuscript/mps_manuscript/figures"

## ---- verified data ----
g65682 <- readRDS(file.path(OUT, "mps_GSE65682_full.rds"))
imm <- readRDS("H:/data/.openclaw/workspace/sepsis_project/output/03_immunosenescence/02_scoring/scores_GSE65682.rds")
d <- merge(g65682[, c("sample_id","status","time","mps","ImmuneSensing","DeathExecution")],
           imm[, c("sample_id","immunity_score","age","sex")], by="sample_id")
d <- d[complete.cases(d), ]
gz <- readRDS(file.path(OUT, "mps_GSE95233_scores.rds")); dz <- gz$dat
dav <- readRDS(file.path(OUT, "mps_Davenport_srs.rds"))
f2 <- readRDS(file.path(OUT, "mps_fig2_data.rds"))
genes2 <- c("MYD88","TLR4","NINJ1","BAK1","BAX")

m1 <- coxph(Surv(time, status) ~ scale(mps), data=d)
m2 <- coxph(Surv(time, status) ~ scale(mps) + age + sex, data=d)
m3 <- coxph(Surv(time, status) ~ scale(mps) + scale(immunity_score) + age + sex, data=d)

## ---- Fig1E: three-cohort forest ----
pdf(file.path(FIG, "Fig1E_crosscohort.pdf"), width=7, height=4)
par(mar=c(4.5,6,3,1))
cc <- data.frame(cohort=c("GSE65682 (n=479)","GSE95233 (n=102)","Davenport (n=265)"),
                 est=c(0.772, 0.449, 1.380), lo=c(0.651, 0.268, 0.994), hi=c(0.915, 0.708, 1.981),
                 metric=c("HR per SD","OR per SD","OR per SD"), p=c(0.003, 0.001, 0.068))
plot(NA, xlim=c(0.2,2.2), ylim=c(0.5,3.5), xlab="HR / OR per SD (28-day mortality)", ylab="", yaxt="n",
     main="Cross-cohort MPS effect")
abline(v=1, lty=2, col="grey50")
for (i in 1:3) {
  yy <- 3.5 - i
  col <- ifelse(cc$est[i]<1, "steelblue", "firebrick")
  segments(cc$lo[i], yy, cc$hi[i], yy, lwd=3, col=col)
  points(cc$est[i], yy, pch=18, cex=2, col=col)
  txt <- sprintf("%s  %s=%.2f (%.2f-%.2f), p=%.3g", cc$cohort[i], cc$metric[i], cc$est[i], cc$lo[i], cc$hi[i], cc$p[i])
  text(cc$lo[i]-0.04, yy, txt, adj=1, cex=0.85)
}
mtext("Davenport divergent direction: endotype-dependent (Fig 5)", side=3, line=0.5, cex=0.8, col="firebrick")
dev.off()

## ---- Fig2: single genes, three cohorts, three points per gene ----
pdf(file.path(FIG, "Fig2_single_genes.pdf"), width=10, height=6.5)
par(mar=c(4.5,6.5,3,1))
plot(NA, xlim=c(0.15,3.2), ylim=c(0.3, length(genes2)+0.8), xlab="Effect size (HR / OR / fold change)", ylab="", yaxt="n",
     main="Single-gene cross-cohort consistency", log="x")
abline(v=1, lty=2, col="grey50")
cols <- c("steelblue","darkgreen","firebrick")
pchs <- c(18, 17, 15)
labs <- c("GSE65682 HR","GSE95233 OR","GSE13904 FC (sepsis vs control)")
for (i in seq_along(genes2)) {
  gn <- genes2[i]
  yy <- length(genes2) - i + 1
  pts <- c(f2$hrs2[gn], f2$or2[gn], f2$fc2[gn])
  for (k in 1:3) {
    if (is.na(pts[k])) next
    points(pts[k], yy, pch=pchs[k], cex=1.7, col=cols[k])
    # annotate value
    vtxt <- sprintf("%.2f", pts[k])
    text(pts[k], yy+0.22, vtxt, cex=0.7, col=cols[k])
  }
  # p-value annotations
  pv <- c(0.0047, 0.275, 0.044, 0.011, 0.229)
  text(0.16, yy, sprintf("%s  (HR p=%.3g | OR p=%.3g | FC p=%.2g)",
       gn, pv[i], c(3.2e-4, 0.885, 0.060, 0.429, 0.121)[i], c(9.5e-7, 9.1e-5, 0.0086, 0.48, 0.45)[i]), adj=0, cex=0.75)
}
legend("bottomright", legend=labs, pch=pchs, col=cols, bty="n", cex=0.9)
dev.off()

## ---- Fig3: three models + C-index ----
pdf(file.path(FIG, "Fig3_immune_independence.pdf"), width=12, height=7)
par(mfrow=c(1,3), mar=c(4.5,4.5,3,1))
# A scatter
plot(d$immunity_score, d$mps, pch=20, col=rgb(0.3,0.4,0.8,0.4), xlab="Pan-immune activation score", ylab="MPS", main="Orthogonality (GSE65682)")
rr <- cor(d$immunity_score, d$mps, method="spearman")
text(par("usr")[1]+0.03, par("usr")[4]-0.06, sprintf("Spearman r = %.3f", rr), adj=0, cex=1.1)
# B three-model forest
models <- list("MPS (unadjusted)"=m1, "+ age, sex"=m2, "+ age, sex, immunity"=m3)
plot(NA, xlim=c(0.5,1.2), ylim=c(0.5,3.5), xlab="MPS HR per SD", ylab="", yaxt="n", main="MPS adjustment sequence")
abline(v=1, lty=2, col="grey50")
for (i in seq_along(models)) {
  mm <- models[[i]]; s <- summary(mm)
  yy <- 3.5 - i
  b <- coef(mm)["scale(mps)"]; ci <- confint(mm)["scale(mps)",]
  segments(exp(ci[1]), yy, exp(ci[2]), yy, lwd=3, col="steelblue")
  points(exp(b), yy, pch=18, cex=2, col="steelblue")
  txt <- sprintf("%s: HR=%.3f (%.3f-%.3f), p=%.4g", names(models)[i], exp(b), exp(ci[1]), exp(ci[2]), s$coefficients["scale(mps)","Pr(>|z|)"])
  text(0.52, yy, txt, adj=0, cex=0.85)
}
# C C-index bars (coxph concordance)
cbase <- summary(coxph(Surv(time, status) ~ age + sex + scale(immunity_score), data=d))$concordance[1]
cnew <- summary(coxph(Surv(time, status) ~ age + sex + scale(immunity_score) + scale(mps), data=d))$concordance[1]
ca <- summary(m1)$concordance[1]
bp <- barplot(c(ca, cbase, cnew), names.arg=c("MPS alone","age+sex+immune","+ MPS"), col=c("steelblue","grey70","steelblue"),
              ylim=c(0,0.9), main="C-index (28-day mortality)", ylab="C-index")
text(bp, c(ca, cbase, cnew)+0.03, sprintf("%.3f", c(ca, cbase, cnew)))
text(2, 0.85, sprintf("delta vs base = +%.3f", cnew-cbase), cex=1.0, font=2)
dev.off()
cat("done\n")
