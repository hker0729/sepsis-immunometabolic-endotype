# Pneumonia-adjusted sensitivity analysis for quadrant Cox models (GSE65682)
# Q1 vs rest and Q3/Q4 vs Q1, with and without pneumonia covariate
suppressMessages({library(survival)})
scores <- readRDS("H:/data/.openclaw/workspace/sepsis_project/output/03_immunosenescence/02_scoring/scores_GSE65682.rds")
cat("cols check: pneumonia_diagnoses:", "pneumonia_diagnoses" %in% names(scores),
    "| quad_short:", "quad_short" %in% names(scores),
    "| mortality:", "mortality_event_28days" %in% names(scores), "\n")
cat("pneumonia NA:", sum(is.na(scores$pneumonia_diagnoses)), "/", nrow(scores), "\n")
cat("pneumonia table:\n"); print(table(scores$pneumonia_diagnoses, useNA="ifany"))

d <- data.frame(
  time = as.numeric(as.character(scores$time_to_event_28days)),
  event = as.numeric(scores$mortality_event_28days),
  age = as.numeric(scores$age),
  sex = as.factor(scores$gender),
  pna_raw = as.character(scores$pneumonia_diagnoses),
  q = scores$quad_short,
  stringsAsFactors = FALSE
)
d$pna_raw[is.na(d$pna_raw) | d$pna_raw == "NA"] <- "none"
d$pna_any <- ifelse(d$pna_raw %in% c("cap","hap"), 1, 0)   # any pneumonia vs none
d$pna3 <- factor(d$pna_raw, levels = c("none","cap","hap")) # 3-level
d <- d[complete.cases(d[, c("time","event","age","sex","q")]), ]
cat("\nanalysis n:", nrow(d), "\n")

# --- Model A: Q1 vs rest ---
d$q1 <- ifelse(d$q == "Q1", 1, 0)
m1 <- coxph(Surv(time, event) ~ q1 + age + sex, data = d)
m2 <- coxph(Surv(time, event) ~ q1 + age + sex + pna_any, data = d)
m2b <- coxph(Surv(time, event) ~ q1 + age + sex + pna3, data = d)
cat("\n=== Q1 vs rest ===\n")
cat("age+sex:        HR", exp(coef(m1)["q1"]), "95%CI", exp(confint(m1)["q1",1]), "-", exp(confint(m1)["q1",2]), "p", summary(m1)$coefficients["q1","Pr(>|z|)"], "\n")
cat("age+sex+pna_any: HR", exp(coef(m2)["q1"]), "95%CI", exp(confint(m2)["q1",1]), "-", exp(confint(m2)["q1",2]), "p", summary(m2)$coefficients["q1","Pr(>|z|)"], "\n")
cat("age+sex+pna3:   HR", exp(coef(m2b)["q1"]), "95%CI", exp(confint(m2b)["q1",1]), "-", exp(confint(m2b)["q1",2]), "p", summary(m2b)$coefficients["q1","Pr(>|z|)"], "\n")
cat("pna_any in m2:  HR", exp(coef(m2)["pna_any"]), "p", summary(m2)$coefficients["pna_any","Pr(>|z|)"], "\n")

# --- Model B: Q3/Q4 vs Q1 (as in fig1d) ---
d$q <- factor(d$q, levels = c("Q1","Q2","Q3","Q4"))
m3 <- coxph(Surv(time, event) ~ q + age + sex, data = d)
m4 <- coxph(Surv(time, event) ~ q + age + sex + pna_any, data = d)
m4b <- coxph(Surv(time, event) ~ q + age + sex + pna3, data = d)
cat("\n=== Q3/Q4 vs Q1 (fig1d style) ===\n")
s3 <- summary(m3); s4 <- summary(m4); s4b <- summary(m4b)
for (lv in c("Q2","Q3","Q4")) {
  nm <- paste0("q", lv)
  cat(lv, "age+sex:     HR", round(exp(coef(m3)[nm]),3), "CI", round(exp(confint(m3)[nm,1]),3), "-", round(exp(confint(m3)[nm,2]),3), "p", round(s3$coefficients[nm,"Pr(>|z|)"],4), "\n")
  cat(lv, "age+sex+pna: HR", round(exp(coef(m4)[nm]),3), "CI", round(exp(confint(m4)[nm,1]),3), "-", round(exp(confint(m4)[nm,2]),3), "p", round(s4$coefficients[nm,"Pr(>|z|)"],4), "\n")
  cat(lv, "age+sex+pna3:HR", round(exp(coef(m4b)[nm]),3), "CI", round(exp(confint(m4b)[nm,1]),3), "-", round(exp(confint(m4b)[nm,2]),3), "p", round(s4b$coefficients[nm,"Pr(>|z|)"],4), "\n")
}

# pneumonia distribution by quadrant (is Q1 enriched for pneumonia?)
tab <- table(d$q, d$pna_any)
cat("\nquadrant x any-pneumonia:\n"); print(tab)
cat("pneumonia rate Q1:", round(prop.table(tab,1)["Q1",2]*100,1), "% vs Q2-Q4:", round(prop.table(tab,1)[c("Q2","Q3","Q4"),2]*100,1), "%\n")
cat("Fisher overall:", fisher.test(tab)$p.value, "\n")
cat("mortality by pna: "); print(round(prop.table(table(d$pna_any, d$event),1)*100,1))
