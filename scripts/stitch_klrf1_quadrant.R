# Stitch analysis: KLRF1 x Quadrant (integration of quadrant paper + KLRF1 paper)
suppressMessages({library(survival)})
set.seed(42)
OUT <- "H:/data/.openclaw/workspace/sepsis_project/output/03_immunosenescence/20_stitch"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

scores <- readRDS("H:/data/.openclaw/workspace/sepsis_project/output/03_immunosenescence/02_scoring/scores_GSE65682.rds")
geo <- readRDS("H:/data/.openclaw/workspace/sepsis_project/data/geo/GSE65682_processed_full.rds")
expr <- geo$exprs

# --- build analysis df ---
d <- data.frame(
  sample_id = rownames(scores),
  klrf1 = as.numeric(expr["KLRF1", rownames(scores)]),
  quad = factor(scores$quadrant, levels = c("Q1 Immune-metabolic failure",
                                             "Q2 Metabolic compensation / immune exhaustion",
                                             "Q3 Immune-preserved / metabolic suppression",
                                             "Q4 Immune-metabolic fitness")),
  quad_short = factor(scores$quad_short, levels = c("Q1","Q2","Q3","Q4")),
  quad_num = scores$quad_num,
  time = as.numeric(scores$time_to_event_28days),
  event = as.numeric(scores$mortality_event_28days),
  age = as.numeric(scores$age),
  sex = scores$gender,
  pneumonia = scores$pneumonia_diagnoses,
  immunity_score = scores$immunity_score,
  nk_function = scores$Immunity_NK_Function,
  nk_effector = scores$Immunity_NK_Effector
)
cat("=== N / events ===\n")
cat(sprintf("N=%d, events=%d (%.1f%%)\n", nrow(d), sum(d$event), 100*mean(d$event)))
cat("KLRF1 expression range: %.3f - %.3f\n", range(d$klrf1))
cat("quad table:\n"); print(table(d$quad_short))

# --- A. KLRF1 by quadrant ---
cat("\n=== A. KLRF1 expression across quadrants ===\n")
print(tapply(d$klrf1, d$quad_short, function(x) round(mean(x), 3)))
kw <- kruskal.test(klrf1 ~ quad_short, data = d)
cat(sprintf("Kruskal-Wallis p = %.3e\n", kw$p.value))
# pairwise Q1 vs others
for (q in c("Q2","Q3","Q4")) {
  w <- wilcox.test(d$klrf1[d$quad_short=="Q1"], d$klrf1[d$quad_short==q])
  cat(sprintf("Q1 vs %s: Wilcoxon p = %.3e (Q1 mean %.3f vs %s mean %.3f)\n",
              q, w$p.value, mean(d$klrf1[d$quad_short=="Q1"]), q, mean(d$klrf1[d$quad_short==q])))
}

# --- B. KLRF1 continuous Cox (overall) ---
cat("\n=== B. KLRF1 continuous Cox (overall) ===\n")
m1 <- coxph(Surv(time, event) ~ klrf1, data = d)
s1 <- summary(m1)
cat(sprintf("Univariable: HR=%.3f [%.3f-%.3f] p=%.3e (per unit)\n",
            exp(coef(m1)), exp(confint(m1))[1], exp(confint(m1))[2], s1$coefficients[5]))
m2 <- coxph(Surv(time, event) ~ klrf1 + age + sex, data = d)
s2 <- summary(m2)
cat(sprintf("Adjusted(age+sex): HR=%.3f [%.3f-%.3f] p=%.3e\n",
            exp(coef(m2)[1]), exp(confint(m2))[1,1], exp(confint(m2))[1,2], s2$coefficients[1,5]))
d$pn <- factor(ifelse(is.na(d$pneumonia), "none", as.character(d$pneumonia)))
m2b <- coxph(Surv(time, event) ~ klrf1 + age + sex + pn, data = d)
s2b <- summary(m2b)
cat(sprintf("Adjusted(age+sex+pneumonia): HR=%.3f [%.3f-%.3f] p=%.3e\n",
            exp(coef(m2b)[1]), exp(confint(m2b))[1,1], exp(confint(m2b))[1,2], s2b$coefficients[1,5]))

# --- C. KLRF1 x Quadrant interaction Cox ---
cat("\n=== C. KLRF1 x Quadrant interaction (Q1 reference) ===\n")
m3 <- coxph(Surv(time, event) ~ klrf1 * quad_short, data = d)
s3 <- summary(m3)
print(round(s3$coefficients[, c("coef", "exp(coef)", "Pr(>|z|)")], 4))
# joint interaction test (LRT vs main-effects model)
m3_main <- coxph(Surv(time, event) ~ klrf1 + quad_short, data = d)
lrt <- anova(m3_main, m3)
cat(sprintf("Interaction LRT: chi2=%.2f df=%d p=%.3e\n", lrt$Chisq[2], lrt$Df[2], lrt$`Pr(>|Chi|)`[2]))

# per-quadrant KLRF1 HR (stratified)
cat("\nPer-quadrant KLRF1 HR (within each quadrant):\n")
per_q <- data.frame()
for (q in c("Q1","Q2","Q3","Q4")) {
  sub <- d[d$quad_short == q, ]
  mq <- coxph(Surv(time, event) ~ klrf1, data = sub)
  ci <- exp(confint(mq))
  per_q <- rbind(per_q, data.frame(quad = q, n = nrow(sub), events = sum(sub$event),
                                   HR = exp(coef(mq)), lo = ci[1], hi = ci[2],
                                   p = summary(mq)$coefficients[5]))
}
print(per_q, row.names = FALSE)

# --- D. Q1-internal KLRF1 high/low (exploratory) ---
cat("\n=== D. Q1-internal KLRF1 median split ===\n")
d$klrf1_med <- ifelse(d$klrf1 > median(d$klrf1), "high", "low")
for (q in c("Q1","Q2","Q3","Q4")) {
  sub <- d[d$quad_short == q, ]
  tab <- table(sub$klrf1_med)
  mk <- survdiff(Surv(time, event) ~ klrf1_med, data = sub)
  ck <- coxph(Surv(time, event) ~ klrf1_med, data = sub)
  ci <- exp(confint(ck))
  cat(sprintf("%s: n_high=%d n_low=%d, logrank p=%.3f, HR(high vs low)=%.2f [%.2f-%.2f] p=%.3f\n",
              q, tab["high"], tab["low"], 1 - pchisq(mk$chisq, 1),
              exp(coef(ck)), ci[1], ci[2], summary(ck)$coefficients[5]))
}

# --- E. correlation with NK/immunity scores ---
cat("\n=== E. KLRF1 vs immune scores (Spearman) ===\n")
for (sc in c("immunity_score", "nk_function", "nk_effector")) {
  ct <- cor.test(d$klrf1, d[[sc]], method = "spearman")
  cat(sprintf("KLRF1 vs %s: rho=%.3f p=%.3e\n", sc, ct$estimate, ct$p.value))
}

# --- save ---
res <- list(
  kw_p = kw$p.value,
  pairwise = sapply(c("Q2","Q3","Q4"), function(q) wilcox.test(d$klrf1[d$quad_short=="Q1"], d$klrf1[d$quad_short==q])$p.value),
  cox_overall = c(HR = exp(coef(m1)), lo = exp(confint(m1))[1], hi = exp(confint(m1))[2], p = s1$coefficients[5]),
  cox_adj = c(HR = exp(coef(m2)[1]), lo = exp(confint(m2))[1,1], hi = exp(confint(m2))[1,2], p = s2$coefficients[1,5]),
  cox_adj_pn = c(HR = exp(coef(m2b)[1]), lo = exp(confint(m2b))[1,1], hi = exp(confint(m2b))[1,2], p = s2b$coefficients[1,5]),
  interaction_lrt_p = lrt$`Pr(>|Chi|)`[2],
  per_quadrant = per_q
)
saveRDS(res, file.path(OUT, "stitch_results.rds"))
write.csv(per_q, file.path(OUT, "per_quadrant_klrf1_hr.csv"), row.names = FALSE)
cat("\nSaved stitch_results.rds + per_quadrant_klrf1_hr.csv\n")

# --- figures ---
png(file.path(OUT, "fig_klrf1_by_quadrant.png"), width = 2200, height = 1600, res = 300)
par(mar = c(6, 5, 3, 1))
cols <- c("#C0392B","#E67E22","#27AE60","#2980B9")
boxplot(klrf1 ~ quad_short, data = d, col = cols, outline = FALSE,
        xlab = "", ylab = "KLRF1 expression", main = "KLRF1 by immune-metabolic quadrant",
        cex.axis = 1.1, cex.lab = 1.2)
stripchart(klrf1 ~ quad_short, data = d, method = "jitter", jitter = 0.15,
           pch = 16, cex = 0.3, col = adjustcolor("grey30", 0.4), vertical = TRUE, add = TRUE)
text(1, max(d$klrf1), sprintf("KW p = %.1e", kw$p.value), cex = 1.1)
dev.off()

# KM in Q1 (high vs low) + other quadrants combined
png(file.path(OUT, "fig_km_klrf1_q1.png"), width = 2200, height = 1600, res = 300)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))
sub1 <- d[d$quad_short == "Q1", ]
fit1 <- survfit(Surv(time, event) ~ klrf1_med, data = sub1)
plot(fit1, col = c("#C0392B", "#2980B9"), lwd = 2, xlab = "Days", ylab = "Survival",
     main = "Q1 (immunometabolic failure)\nKLRF1 high vs low")
legend("topright", c("KLRF1 high", "KLRF1 low"), col = c("#C0392B", "#2980B9"), lwd = 2, bty = "n")
subO <- d[d$quad_short != "Q1", ]
fitO <- survfit(Surv(time, event) ~ klrf1_med, data = subO)
plot(fitO, col = c("#C0392B", "#2980B9"), lwd = 2, xlab = "Days", ylab = "Survival",
     main = "Q2-Q4 combined\nKLRF1 high vs low")
legend("topright", c("KLRF1 high", "KLRF1 low"), col = c("#C0392B", "#2980B9"), lwd = 2, bty = "n")
dev.off()
cat("Figures saved\n")
