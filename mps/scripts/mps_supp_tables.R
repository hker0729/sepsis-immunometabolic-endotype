# Supplementary tables for MPS manuscript
suppressMessages({library(writexl)})
OUT <- "H:/data/.openclaw/workspace/sepsis_project/ltie_manuscript/mps_preanalysis"
SUP <- "H:/data/.openclaw/workspace/sepsis_project/ltie_manuscript/mps_manuscript"

## Table S1: cohort characteristics
imm <- readRDS("H:/data/.openclaw/workspace/sepsis_project/output/03_immunosenescence/02_scoring/scores_GSE65682.rds")
gz <- readRDS(file.path(OUT, "mps_GSE95233_scores.rds"))
dav <- readRDS(file.path(OUT, "mps_Davenport_srs.rds"))
g4 <- readRDS("H:/data/.openclaw/workspace/sepsis_project/data/geo/GSE13904_processed.rds")

tabS1 <- data.frame(
  Cohort = c("GSE65682 (training)", "GSE95233 (validation)", "E-MTAB-4421 (Davenport)", "GSE13904 (case-control)"),
  Disease = c("ICU sepsis (MARS, Netherlands)", "Septic shock (Italy)", "Community-acquired pneumonia (UK)", "Pediatric SIRS/sepsis (USA)"),
  Platform = c("Affymetrix U219", "GE CodeLink UniSet (GPL4204)", "Illumina HumanHT-12 v4", "Affymetrix HG-U133 Plus 2.0"),
  Sample = c("Whole blood", "Whole blood", "Whole blood", "Whole blood"),
  N_total = c(479, 124, 265, 227),
  N_analyzed = c(479, 102, 265, 227),
  Deaths_28d = c(114, 34, 56, NA),
  Mortality_pct = c(23.8, 33.3, 21.1, NA),
  Age_available = c(479, 102, 265, 0),
  Sex_available = c(479, 102, 265, 0),
  MPS_genes_measured = c("16/17 (RPTOR absent)", "17/17", "17/17", "17/17"),
  Outcome = c("28-day mortality (time-to-event)", "28-day survival status", "28-day survival status", "Sepsis vs control (case-control)")
)

## Table S2: 17 MPS genes and modules
tabS2 <- data.frame(
  Gene = c("TLR2","TLR4","TLR7","TLR8","MYD88","NLRP3","TNFRSF1A","TNFRSF1B","NFE2L2",
           "BAX","BAK1","BID","NINJ1","RHOA","RICTOR","MTOR","RPTOR"),
  Module = c(rep("Immune sensing", 9), rep("Death execution", 8)),
  Role_in_mitoxyperilysis = c(
    "TLR2: innate sensing receptor","TLR4: innate sensing receptor (LPS)","TLR7: endosomal RNA sensing",
    "TLR8: endosomal RNA sensing","MYD88: obligate TLR adaptor","NLRP3: inflammasome sensor",
    "TNFRSF1A: TNF receptor 1","TNFRSF1B: TNF receptor 2","NFE2L2: Nrf2 oxidative stress response",
    "BAX: mitochondrial outer membrane permeabilization","BAK1: mitochondrial outer membrane permeabilization",
    "BID: BH3-only activator","NINJ1: plasma membrane rupture executor","RHOA: cytoskeletal regulator (lamellipodia)",
    "RICTOR: mTORC2 component","MTOR: mTORC2 component","RPTOR: mTORC1 component (regulatory)")
)

## Table S3: MR full results
mr <- read.csv(file.path(OUT, "mps_mr_wald.csv"), stringsAsFactors=FALSE)
tabS3 <- mr[, c("gene","SNP","A1","A2","beta_exposure","se_exposure","beta_outcome","se_outcome","OR","ci_lo","ci_hi","p","F","n")]
if (!all(c("SNP","A1","A2","beta_exposure","se_exposure","beta_outcome","se_outcome") %in% colnames(mr))) {
  # fallback: rebuild from known values if csv lacks columns
  tabS3 <- data.frame(
    gene = mr$gene, OR = mr$OR, ci_lo = mr$ci_lo, ci_hi = mr$ci_hi, p = mr$p, F = mr$F, n = mr$n,
    note = "Full SNP-level columns in mps_mr_wald.csv (source of truth)"
  )
}

## Table S4: Davenport SRS interaction complete output
dav$srs2f <- factor(dav$srs2, levels=c(0,1), labels=c("SRS1","SRS2"))
m_full <- glm(dead ~ scale(MPS)*srs2f, data=dav, family=binomial)
s_full <- summary(m_full)
tabS4 <- data.frame(
  Term = rownames(s_full$coefficients),
  Beta = s_full$coefficients[,1],
  SE = s_full$coefficients[,2],
  OR = exp(s_full$coefficients[,1]),
  p = s_full$coefficients[,4]
)
tabS4b <- data.frame(
  Stratum = c("SRS1","SRS2","Interaction"),
  N = c(108, 157, 265),
  Deaths = c(29, 27, 56),
  MPS_OR_per_SD = c(1.035, 1.570, NA),
  CI95 = c("0.677-1.616", "0.964-2.750", NA),
  p = c(0.875, 0.093, 0.308)
)

## Table S5: single-gene cross-cohort effects
f2 <- readRDS(file.path(OUT, "mps_fig2_data.rds"))
tabS5 <- data.frame(
  Gene = c("MYD88","TLR4","NINJ1","BAK1","BAX"),
  GSE65682_HR = unname(f2$hrs2), GSE65682_p = c(0.0047, 0.275, 0.044, 0.011, 0.229),
  GSE95233_OR = unname(f2$or2), GSE95233_p = c(3.2e-4, 0.885, 0.060, 0.429, 0.121),
  GSE13904_FC = unname(f2$fc2), GSE13904_p = c(9.5e-7, 9.1e-5, 0.0086, 0.48, 0.45)
)

write_xlsx(list(
  "TableS1_Cohorts" = tabS1,
  "TableS2_MPS_genes" = tabS2,
  "TableS3_MR_results" = tabS3,
  "TableS4_SRS_interaction" = tabS4b,
  "TableS5_Single_gene" = tabS5
), file.path(SUP, "Supplementary_Tables.xlsx"))
cat("Supplementary_Tables.xlsx written\n")
