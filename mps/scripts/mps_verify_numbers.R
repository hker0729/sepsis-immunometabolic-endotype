# Verify all numbers needed for manuscript (anti-fabrication check)
OUT <- "H:/data/.openclaw/workspace/sepsis_project/ltie_manuscript/mps_preanalysis"

cat("===== 1. MR Wald results =====\n")
mr <- read.csv(file.path(OUT, "mps_mr_wald.csv"), stringsAsFactors=FALSE)
print(mr)

cat("\n===== 2. MYD88 single-gene Cox in GSE65682 =====\n")
suppressMessages(library(survival))
mps <- readRDS(file.path(OUT, "mps_GSE65682_scores.rds"))
cat("cols:", paste(colnames(mps), collapse=" | "), "\n")
cat("n =", nrow(mps), "\n")
# single gene cox for MYD88, TLR4, NINJ1, BAK1, BAX
for (g in c("MYD88","TLR4","NINJ1","BAK1","BAX")) {
  if (!g %in% colnames(mps)) { cat("  ", g, "not in cols\n"); next }
  m <- coxph(Surv(time, status) ~ scale(mps[[g]]), data=mps)
  s <- summary(m)
  cat(sprintf("  %-8s HR=%.3f (%.3f-%.3f) p=%.4g\n", g, exp(coef(m)), exp(confint(m)[,1]), exp(confint(m)[,2]), s$coefficients[,"Pr(>|z|)"]))
}
# module correlation
cat("\n===== 3. ImmuneSensing vs DeathExecution correlation (GSE65682) =====\n")
if (all(c("ImmuneSensing","DeathExecution") %in% colnames(mps))) {
  r <- cor(mps$ImmuneSensing, mps$DeathExecution, method="spearman")
  cat(sprintf("  r = %.3f\n", r))
} else cat("  module cols not in mps_GSE65682_scores.rds\n")

cat("\n===== 4. Module Cox (GSE65682) =====\n")
for (mod in c("ImmuneSensing","DeathExecution")) {
  if (!mod %in% colnames(mps)) next
  m <- coxph(Surv(time, status) ~ scale(mps[[mod]]), data=mps)
  s <- summary(m)
  cat(sprintf("  %-15s HR=%.3f (%.3f-%.3f) p=%.4g\n", mod, exp(coef(m)), exp(confint(m)[,1]), exp(confint(m)[,2]), s$coefficients[,"Pr(>|z|)"]))
}

cat("\n===== 5. GSE13904 MYD88 FC =====\n")
g <- readRDS("H:/data/.openclaw/workspace/sepsis_project/data/geo/GSE13904_processed.rds")
ex <- g$exprs; pd <- g$pdata
my <- ex["MYD88", ]; tl <- ex["TLR4", ]
cat(sprintf("  MYD88 sepsis mean=%.3f control mean=%.3f FC=%.2f\n", mean(my[pd$group=="sepsis"]), mean(my[pd$group=="control"]), mean(my[pd$group=="sepsis"])/mean(my[pd$group=="control"])))
cat(sprintf("  TLR4  sepsis mean=%.3f control mean=%.3f FC=%.2f\n", mean(tl[pd$group=="sepsis"]), mean(tl[pd$group=="control"]), mean(tl[pd$group=="sepsis"])/mean(tl[pd$group=="control"])))

cat("\n===== 6. GSE65682 clinical vars for IDI/NRI =====\n")
imm <- readRDS("H:/data/.openclaw/workspace/sepsis_project/output/03_immunosenescence/02_scoring/scores_GSE65682.rds")
cat("imm cols:", paste(colnames(imm), collapse=" | "), "\n")
cat("age NA:", sum(is.na(imm$age)), "| sex NA:", sum(is.na(imm$sex)), "\n")
cat("mps cols check - has immune module?\n")
