# Fig 1d: multivariable Cox forest (Q1 as reference + Q1 vs Q2-4), adjusted age+sex
suppressMessages(library(survival))
scores <- readRDS("H:/data/.openclaw/workspace/sepsis_project/output/03_immunosenescence/02_scoring/scores_GSE65682.rds")
d <- data.frame(
  time = as.numeric(scores$time_to_event_28days),
  event = as.numeric(scores$mortality_event_28days),
  age = as.numeric(scores$age),
  sex = scores$gender,
  quad = factor(scores$quad_short, levels = c("Q1","Q2","Q3","Q4"))
)
d$q1 <- ifelse(d$quad == "Q1", 1, 0)

# model A: Q1 vs Q2-4
mA <- coxph(Surv(time, event) ~ q1 + age + sex, data = d)
ciA <- exp(confint(mA))[1, ]
# model B: quadrant with Q1 ref
mB <- coxph(Surv(time, event) ~ quad + age + sex, data = d)
ciB <- exp(confint(mB))[1:3, ]

rows <- rbind(
  c("Q1 vs Q2-Q4 (adj)", exp(coef(mA)[1]), ciA[1], ciA[2], summary(mA)$coefficients[1,5]),
  c("Q2 vs Q1", exp(coef(mB)[1]), ciB[1,1], ciB[1,2], summary(mB)$coefficients[1,5]),
  c("Q3 vs Q1", exp(coef(mB)[2]), ciB[2,1], ciB[2,2], summary(mB)$coefficients[2,5]),
  c("Q4 vs Q1", exp(coef(mB)[3]), ciB[3,1], ciB[3,2], summary(mB)$coefficients[3,5])
)
colnames(rows) <- c("term","HR","lo","hi","p")
rows <- as.data.frame(rows, stringsAsFactors = FALSE)
for (i in 2:5) rows[[i]] <- as.numeric(rows[[i]])
rows$p <- format.pval(rows$p, digits = 3, eps = 1e-4)
print(rows, row.names = FALSE)

OUT <- "H:/data/.openclaw/workspace/sepsis_project/output/03_immunosenescence/20_stitch"
png(file.path(OUT, "fig1d_quadrant_adjusted_forest.png"), width = 2200, height = 1400, res = 300)
par(mar = c(5, 8, 3, 2))
n <- nrow(rows)
x <- log(rows$HR); lo <- log(rows$lo); hi <- log(rows$hi)
plot(x, n:1, xlim = log(c(0.2, 3.2)), ylim = c(0.5, n + 0.5), pch = 15, cex = 1.4,
     col = c("#C0392B", rep("#34495E", 3)), xaxt = "n", yaxt = "n",
     xlab = "Hazard ratio (95% CI, log scale)", ylab = "", main = "28-day mortality: multivariable Cox (adjusted age + sex)")
abline(v = log(1), lty = 2, col = "grey50")
for (i in 1:n) {
  lines(c(lo[i], hi[i]), rep(n - i + 1, 2), lwd = 2.5, col = c("#C0392B", rep("#34495E", 3))[i])
}
axis(1, at = log(c(0.25, 0.5, 1, 2, 3)), labels = c(0.25, 0.5, 1, 2, 3))
text(log(0.2), n:1 + 0.25, adj = 0, cex = 0.85,
     labels = sprintf("%s: HR=%.2f [%.2f-%.2f], p=%s", rows$term, rows$HR, rows$lo, rows$hi, rows$p))
dev.off()
write.csv(rows, file.path(OUT, "fig1d_forest_data.csv"), row.names = FALSE)
cat("\nForest saved.\n")
