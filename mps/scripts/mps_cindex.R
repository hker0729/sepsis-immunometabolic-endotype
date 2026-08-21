# C-index increment (IDI/NRI alternative) + single-cell object info
suppressMessages({library(survival)})
OUT <- "H:/data/.openclaw/workspace/sepsis_project/ltie_manuscript/mps_preanalysis"
pd <- readRDS(file.path(OUT, "mps_GSE65682_full.rds"))
imm <- readRDS("H:/data/.openclaw/workspace/sepsis_project/output/03_immunosenescence/02_scoring/scores_GSE65682.rds")
d <- merge(pd[, c("sample_id","status","time","mps")],
           imm[, c("sample_id","immunity_score","age","sex")], by="sample_id")
d <- d[complete.cases(d), ]
cat("n =", nrow(d), "\n")
cbase <- concordance(Surv(time, status) ~ age + sex + immunity_score, data=d)
cnew  <- concordance(Surv(time, status) ~ age + sex + immunity_score + mps, data=d)
cat(sprintf("C-index base (age+sex+immunity) = %.4f\n", cbase$concordance))
cat(sprintf("C-index +MPS = %.4f  (delta %.4f)\n", cnew$concordance, cnew$concordance - cbase$concordance))
# also MPS alone
cmps <- concordance(Surv(time, status) ~ mps, data=d)
cat(sprintf("C-index MPS alone = %.4f\n", cmps$concordance))

cat("\n===== single-cell object check =====\n")
sc <- readRDS("H:/data/.openclaw/workspace/sepsis_project/output/08_SCISSORS_lilrb2/monocyte_seurat.rds")
cat("monocyte_seurat: cells =", ncol(sc), "| features =", nrow(sc), "\n")
cat("meta cols:", paste(colnames(sc@meta.data), collapse=" | "), "\n")
if ("group" %in% colnames(sc@meta.data)) print(table(sc@meta.data$group))
if ("celltype" %in% colnames(sc@meta.data)) print(table(sc@meta.data$celltype))
# full dataset?
full <- "H:/data/.openclaw/workspace/sepsis_project/output/06_scRNA/seurat_integrated.rds"
if (file.exists(full)) {
  f <- readRDS(full)
  cat("\nfull seurat: cells =", ncol(f), "\n")
  print(table(f@meta.data$celltype))
} else cat("\nno full seurat at", full, "\n")
