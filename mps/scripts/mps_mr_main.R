# MPS MR: NINJ1/BAK1/BAX top cis-eQTL x ieu-b-4980 (Wald ratio, ltie_12 pipeline)
# Requires valid IEU JWT (OPENGWAS_JWT); token refresh needed 2026-08-21
tok <- readLines("H:/data/.openclaw/workspace/sepsis_project/references/ieu_opengwas_jwt.txt", warn=FALSE)
tok <- trimws(paste(tok, collapse=""))
Sys.setenv(OPENGWAS_JWT = tok)
options(ieugwasr.api = "https://api.opengwas.io")
options(ieugwasr.token = tok)
suppressMessages(library(ieugwasr))

OUT <- "H:/data/.openclaw/workspace/sepsis_project/ltie_manuscript/mps_preanalysis"
eq <- read.csv(file.path(OUT, "mps_eqtlgen_subset.csv"), check.names=FALSE, stringsAsFactors=FALSE)
eq$SNP <- gsub('"', '', eq$SNP); eq$AssessedAllele <- gsub('"', '', eq$AssessedAllele)
eq$OtherAllele <- gsub('"', '', eq$OtherAllele); eq$GeneSymbol <- gsub('"', '', eq$GeneSymbol)
eq$Zscore <- as.numeric(eq$Zscore); eq$NrSamples <- as.numeric(eq$NrSamples)

genes <- c("NINJ1","BAK1","BAX")
ivs <- do.call(rbind, lapply(genes, function(g) {
  d <- eq[eq$GeneSymbol == g, ]
  d[order(-abs(d$Zscore)), ][1, ]
}))
cat("top instruments:\n")
print(ivs[, c("GeneSymbol","SNP","AssessedAllele","OtherAllele","Zscore","NrSamples")])

res <- associations(variants = ivs$SNP, id = "ieu-b-4980", proxies = 0)
res$rsid <- as.character(res$rsid)
ivs$SNP <- as.character(ivs$SNP)
m <- merge(ivs, res, by.x="SNP", by.y="rsid", all.x=TRUE)
m <- m[!is.na(m$beta), ]
cat("\nmatched:", nrow(m), "of", nrow(ivs), "\n")

out <- do.call(rbind, lapply(seq_len(nrow(m)), function(i) {
  r <- m[i, ]
  f <- r$eaf; if (is.na(f) || f <= 0 || f >= 1) f <- 0.5
  b_e <- r$Zscore / sqrt(2 * f * (1 - f) * (r$NrSamples + r$Zscore^2))
  if (r$AssessedAllele == r$ea) s <- 1
  else if (r$AssessedAllele == r$nea) s <- -1
  else s <- NA
  if (is.na(s)) return(NULL)
  b_o <- s * r$beta
  b_r <- b_o / b_e; se_r <- r$se / abs(b_e)
  z <- b_r / se_r; p <- 2 * pnorm(-abs(z))
  fstat <- b_e^2 / (r$se / abs(b_e))^2  # approx F
  data.frame(gene = r$GeneSymbol, snp_eqtl = r$SNP, snp_out = r$rsid,
             a1 = r$AssessedAllele, a2 = r$OtherAllele, ea_out = r$ea,
             beta_eqtl = round(b_e, 5), beta_out = round(b_o, 5),
             OR = round(exp(b_r), 3), ci_lo = round(exp(b_r - 1.96*se_r), 3),
             ci_hi = round(exp(b_r + 1.96*se_r), 3),
             p = formatC(p, format="e", digits=2), F = round(fstat, 1),
             n = r$n, stringsAsFactors=FALSE)
}))
cat("\n========== WALD RATIO MR: 表达 -> 脓毒症 (ieu-b-4980) ==========\n")
print(out, row.names=FALSE)
write.csv(out, file.path(OUT, "mps_mr_wald.csv"), row.names=FALSE)

cat("\n方向判定（OR>1 = 表达↑→脓毒症风险↑）:\n")
for (i in seq_len(nrow(out))) {
  o <- out[i, ]
  cat(sprintf("  %s: OR=%.2f p=%s F=%.1f %s\n", o$gene, o$OR, o$p, o$F,
              ifelse(as.numeric(o$p) < 0.05, "**显著**",
                     ifelse(o$OR > 1, "风险方向", "保护方向"))))
}
