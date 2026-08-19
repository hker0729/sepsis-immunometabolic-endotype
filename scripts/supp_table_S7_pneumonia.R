suppressMessages(library(survival))
scores <- readRDS("H:/data/.openclaw/workspace/sepsis_project/output/03_immunosenescence/02_scoring/scores_GSE65682.rds")
pna <- as.character(scores$pneumonia_diagnoses)
pna[is.na(pna) | pna == "NA"] <- "none"
quad <- scores$quad_short
tab3 <- table(quad, factor(pna, levels = c("none", "cap", "hap")))
anyp <- ifelse(pna %in% c("cap", "hap"), "any", "none")
tab2 <- table(quad, anyp)
cat("=== any-pneumonia counts (2x4) ===\n"); print(tab2)
cat("Fisher(any) p:", fisher.test(tab2)$p.value, "\n")
cat("any-pneumonia % by quadrant:", round(100 * prop.table(tab2, 1)[, "any"], 1), "\n")
cat("quadrant n:", table(quad), "\n")
mort <- as.numeric(scores$mortality_event_28days)
cat("death rate no-pna vs pna:", round(100 * tapply(mort, anyp, mean, na.rm = TRUE), 1), "\n")
out <- data.frame(
  quadrant = rownames(tab3),
  n = rowSums(tab3),
  none_n = tab3[, "none"], cap_n = tab3[, "cap"], hap_n = tab3[, "hap"],
  any_pna_n = tab3[, "cap"] + tab3[, "hap"],
  any_pna_pct = round(100 * (tab3[, "cap"] + tab3[, "hap"]) / rowSums(tab3), 1),
  deaths = as.integer(tapply(mort, quad, sum, na.rm = TRUE)),
  mortality_pct = round(100 * tapply(mort, quad, mean, na.rm = TRUE), 1)
)
write.csv(out, "supplementary/S7_pneumonia_by_quadrant.csv", row.names = FALSE)
cat("csv saved\n")
