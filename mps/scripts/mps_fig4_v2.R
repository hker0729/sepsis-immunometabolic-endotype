# Fig4 v2: 4 panels matching revised legend
# A: NINJ1 by group (monocytes) | B: MPS myeloid vs lymphoid (21k PBMC)
# C: Sepsis T0 vs T6 | D: Septic shock T0 vs T6
suppressMessages({library(Seurat)})
FIG <- "H:/data/.openclaw/workspace/sepsis_project/ltie_manuscript/mps_manuscript/figures"
sc <- readRDS("H:/data/.openclaw/workspace/sepsis_project/output/08_SCISSORS_lilrb2/monocyte_seurat.rds")
md <- sc@meta.data
ex <- GetAssayData(sc, assay="RNA", layer="data")
md$NINJ1 <- if ("NINJ1" %in% rownames(ex)) as.numeric(ex["NINJ1", ]) else NA

s21 <- readRDS("H:/data/.openclaw/workspace/sepsis_project/ltie_manuscript/mps_preanalysis/mps_sc_GSE167363.rds")
md21 <- s21@meta.data
myeloid <- c("CD14_Mono","Neutrophil"); lymphoid <- c("NK_cell","CD8_T_cell","CD4_T_cell","B_cell")
my <- md21$celltype %in% myeloid; ly <- md21$celltype %in% lymphoid
w_my <- wilcox.test(md21$MPS[my], md21$MPS[ly])
cat(sprintf("MPS myeloid n=%d mean=%.3f vs lymphoid n=%d mean=%.3f, p=%.3g\n",
    sum(my), mean(md21$MPS[my]), sum(ly), mean(md21$MPS[ly]), w_my$p.value))

sub <- md[md$group=="Sepsis" & !is.na(md$time), ]
sub2 <- md[md$group=="SS" & !is.na(md$time), ]
w1 <- wilcox.test(sub$NINJ1[sub$time=="T0"], sub$NINJ1[sub$time=="T6"])
w2 <- wilcox.test(sub2$NINJ1[sub2$time=="T0"], sub2$NINJ1[sub2$time=="T6"])
cat(sprintf("Sepsis NINJ1 T0 %.3f (n=%d) vs T6 %.3f (n=%d), p=%.4g\n",
    mean(sub$NINJ1[sub$time=="T0"]), sum(sub$time=="T0"), mean(sub$NINJ1[sub$time=="T6"]), sum(sub$time=="T6"), w1$p.value))
cat(sprintf("SS NINJ1 T0 %.3f vs T6 %.3f, p=%.4g\n", mean(sub2$NINJ1[sub2$time=="T0"]), mean(sub2$NINJ1[sub2$time=="T6"]), w2$p.value))

pdf(file.path(FIG, "Fig4_NINJ1_temporal.pdf"), width=11, height=7)
par(mfrow=c(1,4), mar=c(4.5,4.5,3,1))
# A: NINJ1 by group
cols <- c("Healthy"="grey60","NSES"="grey40","Sepsis"="steelblue","SS"="firebrick")
boxplot(NINJ1 ~ group, data=md, col=cols[levels(factor(md$group))], ylab="NINJ1 log-normalized", main="Monocyte NINJ1 by group", outline=FALSE)
# B: MPS myeloid vs lymphoid
boxplot(MPS ~ celltype, data=md21[md21$celltype %in% c(myeloid, lymphoid), ], col=c("steelblue","steelblue","grey70","grey70","grey70","grey70")[match(levels(factor(md21$celltype[md21$celltype %in% c(myeloid,lymphoid)])), c("CD14_Mono","Neutrophil","NK_cell","CD8_T_cell","CD4_T_cell","B_cell"))],
        ylab="MPS", main=sprintf("MPS by cell type\nmyeloid vs lymphoid p=%.1g", w_my$p.value), outline=FALSE, las=2)
# C: Sepsis T0 vs T6
boxplot(NINJ1 ~ time, data=sub, col=c("steelblue","lightblue"), ylab="NINJ1", main=sprintf("Sepsis: T0 vs T6\np=%.4g", w1$p.value))
# D: SS T0 vs T6
boxplot(NINJ1 ~ time, data=sub2, col=c("firebrick","salmon"), ylab="NINJ1", main=sprintf("Septic shock: T0 vs T6\np=%.4g", w2$p.value))
dev.off()
cat("Fig4 v2 written\n")
