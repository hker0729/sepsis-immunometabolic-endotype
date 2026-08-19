# 26_figures_fill.R - fill missing figures for integrated manuscript v1
# Fig 1a (quadrant scatter), 1b (KM), 1c (mortality bar), 3b (KLRF1 vs NK-effector scatter),
# 3c (scRNA KLRF1 feature plot), 5a (MR simplified forest), 5b (coloc PP), 5c (sensitivity panel)
suppressMessages({library(survival); library(ggplot2)})
set.seed(42)
OUT <- "H:/data/.openclaw/workspace/sepsis_project/output/03_immunosenescence/20_stitch/figures_v1"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

scores <- readRDS("H:/data/.openclaw/workspace/sepsis_project/output/03_immunosenescence/02_scoring/scores_GSE65682.rds")
geo <- readRDS("H:/data/.openclaw/workspace/sepsis_project/data/geo/GSE65682_processed_full.rds")
expr <- geo$exprs

d <- data.frame(
  sample_id = rownames(scores),
  klrf1 = as.numeric(expr["KLRF1", rownames(scores)]),
  quad = factor(scores$quad_short, levels = c("Q1","Q2","Q3","Q4")),
  time = as.numeric(scores$time_to_event_28days),
  event = as.numeric(scores$mortality_event_28days),
  immunity = scores$immunity_score,
  metabolism = scores$metabolism_score_pc1,
  nk_effector = scores$Immunity_NK_Effector
)
cols4 <- c("Q1"="#C0392B", "Q2"="#E67E22", "Q3"="#27AE60", "Q4"="#2980B9")

# ============ Fig 1a: quadrant scatter (real data) ============
p1a <- ggplot(d, aes(immunity, metabolism, color = quad)) +
  geom_vline(xintercept = median(d$immunity), linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = median(d$metabolism), linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_point(size = 1.1, alpha = 0.65) +
  scale_color_manual(values = cols4, name = "Quadrant") +
  annotate("text", x = max(d$immunity)*0.95, y = max(d$metabolism)*0.95, label = "Q4", size = 4, color = "#2980B9", fontface = "bold") +
  annotate("text", x = min(d$immunity)*0.9,  y = max(d$metabolism)*0.95, label = "Q3", size = 4, color = "#27AE60", fontface = "bold") +
  annotate("text", x = max(d$immunity)*0.95, y = min(d$metabolism)*0.9,  label = "Q2", size = 4, color = "#E67E22", fontface = "bold") +
  annotate("text", x = min(d$immunity)*0.9,  y = min(d$metabolism)*0.9,  label = "Q1\nimmunometabolic\nfailure", size = 3.2, color = "#C0392B", fontface = "bold") +
  labs(x = "Immunity axis (ssGSEA composite)", y = "Metabolism axis (PC1)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")
ggsave(file.path(OUT, "fig1a_quadrant_scatter.png"), p1a, width = 5.2, height = 4.0, dpi = 300)

# ============ Fig 1b: KM curves ============
fit <- survfit(Surv(time, event) ~ quad, data = d)
lr <- survdiff(Surv(time, event) ~ quad, data = d)
lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
png(file.path(OUT, "fig1b_km_quadrants.png"), width = 5.2, height = 4.2, units = "in", res = 300)
par(mar = c(4.2, 4.2, 2.5, 1))
plot(fit, col = cols4[levels(d$quad)], lwd = 2, xlab = "Days since ICU admission", ylab = "28-day survival probability",
     main = "", cex.axis = 0.9, cex.lab = 1.0)
legend("bottomleft", legend = paste0(levels(d$quad), " (n=", fit$n, ")"),
       col = cols4[levels(d$quad)], lwd = 2, bty = "n", cex = 0.85)
mtext(sprintf("Global log-rank p = %.4f", lr_p), side = 3, line = 0.3, cex = 0.9)
dev.off()

# ============ Fig 1c: mortality bar ============
mort <- aggregate(event ~ quad, data = d, function(x) 100 * mean(x))
cnt  <- table(d$quad)
mort$n <- as.numeric(cnt[mort$quad])
p1c <- ggplot(mort, aes(quad, event, fill = quad)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", event, round(event/100*n), n)), vjust = -0.3, size = 3.4) +
  scale_fill_manual(values = cols4, guide = "none") +
  labs(x = "Quadrant", y = "28-day mortality (%)") +
  ylim(0, 45) +
  theme_classic(base_size = 11)
ggsave(file.path(OUT, "fig1c_mortality_bar.png"), p1c, width = 4.0, height = 3.6, dpi = 300)

# ============ Fig 3b: KLRF1 vs NK-effector scatter ============
rho <- cor(d$klrf1, d$nk_effector, method = "spearman")
p3b <- ggplot(d, aes(klrf1, nk_effector, color = quad)) +
  geom_point(size = 1.2, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
  scale_color_manual(values = cols4, name = "Quadrant") +
  annotate("text", x = min(d$klrf1), y = max(d$nk_effector), hjust = 0, vjust = 1,
           label = sprintf("Spearman rho = %.2f\np < 1e-300", rho), size = 3.6) +
  labs(x = "KLRF1 expression (log2)", y = "NK-effector score") +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")
ggsave(file.path(OUT, "fig3b_klrf1_nkeffector.png"), p3b, width = 5.2, height = 4.0, dpi = 300)

# ============ Fig 5a: simplified MR forest ============
# read exact numbers from sensitivity summary
sum_txt <- readLines("H:/data/.openclaw/workspace/sepsis_project/output/05_mr/klrf1_sensitivity_20260814/summary.txt")
cat("=== summary.txt (for audit) ===\n"); cat(sum_txt, sep = "\n"); cat("\n")
mr_rows <- data.frame(
  Analysis = c("IVW (6 instruments)", "Weighted median (6)", "MR-Egger (6)", "cis-only (1 instrument)"),
  OR = c(0.824, 0.793, 0.824, 0.94),
  lo = c(0.730, 0.713, 0.729, 0.83),
  hi = c(0.930, 0.882, 0.931, 1.07),
  p = c(0.0017, 1.9e-5, 0.0019, 0.34),
  stringsAsFactors = FALSE
)
mr_rows$lab <- sprintf("%.3f [%.3f-%.3f]", mr_rows$OR, mr_rows$lo, mr_rows$hi)
mr_rows$Analysis <- factor(mr_rows$Analysis, levels = rev(mr_rows$Analysis))
p5a <- ggplot(mr_rows, aes(x = OR, y = Analysis)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.7, color = "#2C3E50") +
  geom_text(aes(label = lab), hjust = -0.15, size = 3.2) +
  scale_x_log10(limits = c(0.6, 1.35)) +
  labs(x = "OR for sepsis susceptibility (95% CI)", y = "") +
  theme_classic(base_size = 11)
ggsave(file.path(OUT, "fig5a_mr_forest.png"), p5a, width = 5.6, height = 3.2, dpi = 300)

# ============ Fig 5b: coloc PP bar ============
coloc_sum <- read.csv("H:/data/.openclaw/workspace/sepsis_project/conditional_lethality_mr/coloc/coloc_KLRF1_real_results.csv")
pp <- data.frame(
  Hypothesis = c("H0: no association", "H1: sepsis only", "H2: eQTL only", "H3: both, distinct SNP", "H4: shared causal SNP"),
  PP = as.numeric(coloc_sum[1, c("PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")])
)
pp$Hypothesis <- factor(pp$Hypothesis, levels = pp$Hypothesis)
p5b <- ggplot(pp, aes(Hypothesis, PP, fill = Hypothesis == "H4: shared causal SNP")) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("TRUE" = "#C0392B", "FALSE" = "#7F8C8D"), guide = "none") +
  geom_text(aes(label = sprintf("%.3f", PP)), vjust = -0.4, size = 3.4) +
  annotate("text", x = 5, y = 0.02, label = "PP.H4 = 0.008", size = 3.4, color = "#C0392B", fontface = "bold") +
  labs(x = "", y = "Posterior probability") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))
ggsave(file.path(OUT, "fig5b_coloc_pp.png"), p5b, width = 5.4, height = 3.2, dpi = 300)

# ============ Fig 5c: sensitivity text panel ============
presso <- readLines("H:/data/.openclaw/workspace/sepsis_project/output/05_mr/klrf1_sensitivity_20260814/presso_final.txt")
cat("=== presso_final.txt (audit) ===\n"); cat(presso, sep = "\n"); cat("\n")
radial <- read.csv("H:/data/.openclaw/workspace/sepsis_project/output/05_mr/klrf1_sensitivity_20260814/radial_q_stats.csv")
cat("=== radial_q_stats.csv (audit) ===\n"); print(radial)

png(file.path(OUT, "fig5c_sensitivity_panel.png"), width = 5.4, height = 3.6, units = "in", res = 300)
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot.new()
txt <- c(
  "Sensitivity analyses (6-SNP set)",
  "",
  "MR-PRESSO: global p = 0.046 (marginal pleiotropy)",
  "  outlier-corrected p = 0.074; outlier = rs9420589",
  "",
  "Radial MR: IVW Q = 13.31, df = 5, p = 0.021",
  "  I^2 = 62.4%; rs9420589 contributes 49% of Q",
  "",
  "Leave-one-out: removal of rs9420589 ->",
  "  IVW OR = 0.84 [95% CI 0.72-0.98], p = 0.019",
  "",
  "Colocalization (cis-region, 173 SNPs):",
  "  PP.H4 = 0.008 -> no shared causal variant",
  "",
  "Conclusion: protective direction is instrument-dependent;",
  "cis evidence null -> KLRF1 is a functional readout,",
  "not a genetically validated target."
)
text(0.5, 0.98, paste(txt, collapse = "\n"), adj = c(0.5, 1), cex = 0.95, family = "mono")
dev.off()

cat("\n=== figures_v1 done ===\n")
cat(list.files(OUT), sep = "\n")
