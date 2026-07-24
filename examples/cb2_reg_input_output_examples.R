#!/usr/bin/env Rscript

section <- function(title) {
  cat("\n", strrep("=", 72), "\n", title, "\n",
      strrep("=", 72), "\n", sep = "")
}

if (file.exists(file.path("CB2", "DESCRIPTION")) &&
    requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all("CB2", quiet = TRUE)
} else if (!requireNamespace("CB2", quietly = TRUE)) {
  stop("Install CB2, or run this script from the CB2-Reg repository root.")
}

section("Example 1 input: two experimental groups")

two_group <- data.frame(
  sample = paste0(rep(c("A", "B"), each = 4), 1:4),
  group = factor(rep(c("A", "B"), each = 4)),
  count = c(74, 112, 91, 139, 128, 177, 164, 211),
  total = c(52000, 81000, 69000, 97000, 55000, 76000, 72000, 93000)
)
two_group$proportion <- two_group$count / two_group$total
print(two_group, digits = 9, row.names = FALSE)

fit_a <- fit_ab(
  matrix(two_group$count[two_group$group == "A"], nrow = 1),
  matrix(two_group$total[two_group$group == "A"], nrow = 1)
)
fit_b <- fit_ab(
  matrix(two_group$count[two_group$group == "B"], nrow = 1),
  matrix(two_group$total[two_group$group == "B"], nrow = 1)
)
phat <- c(A = fit_a$phat, B = fit_b$phat)
vhat <- c(A = fit_a$vhat, B = fit_b$vhat)
legacy_t <- diff(phat) / sqrt(sum(vhat))
legacy_df <- sum(vhat)^2 /
  (vhat[1]^2 / 3 + vhat[2]^2 / 3)

section("Example 1A output: original CB2")
original_output <- data.frame(
  phat_A = phat[["A"]],
  phat_B = phat[["B"]],
  effect_B_minus_A = diff(phat),
  std_error = sqrt(sum(vhat)),
  t_value = legacy_t,
  df = legacy_df,
  p_value = 2 * pt(-abs(legacy_t), df = legacy_df)
)
print(original_output, digits = 12, row.names = FALSE)

group_fit <- bbreg(
  count = two_group$count,
  total = two_group$total,
  formula = ~ group,
  data = two_group
)

section("Example 1B output: default CB2-Reg logit model")
print(group_fit$coefficient_table, digits = 12)
cat(
  "The groupB effect is in log-odds units, not raw-proportion units.\n"
)

section("Example 2 input: continuous dose")

dose_input <- data.frame(
  sample = sprintf("dose_%02d", 1:10),
  dose = rep(0:4, each = 2),
  total = rep(c(80000, 100000), 5)
)
dose_input$true_probability <- plogis(-7 + 0.35 * dose_input$dose)
dose_input$count <- round(
  dose_input$total * dose_input$true_probability
)
dose_input$observed_proportion <-
  dose_input$count / dose_input$total
print(dose_input, digits = 9, row.names = FALSE)

dose_fit <- bbreg(
  count = dose_input$count,
  total = dose_input$total,
  formula = ~ dose,
  data = dose_input
)

section("Example 2 output: guide-abundance trend per dose step")
print(dose_fit$coefficient_table, digits = 12)
cat(
  "True dose slope: 0.35; fitted dose slope:",
  sprintf("%.9f", coef(dose_fit)[["dose"]]), "\n"
)

section("Example 3 input: dose with batch adjustment")

adjusted_input <- data.frame(
  sample = sprintf("s%02d", 1:12),
  dose = rep(0:5, 2),
  batch = factor(rep(c("A", "B"), each = 6)),
  total = rep(c(80000, 95000, 85000, 105000, 90000, 100000), 2)
)
adjusted_eta <- -7 + 0.3 * adjusted_input$dose +
  0.2 * (adjusted_input$batch == "B")
adjusted_input$count <- round(
  adjusted_input$total * plogis(adjusted_eta)
)
adjusted_input$observed_proportion <-
  adjusted_input$count / adjusted_input$total
print(adjusted_input, digits = 9, row.names = FALSE)

adjusted_fit <- bbreg(
  count = adjusted_input$count,
  total = adjusted_input$total,
  formula = ~ dose + batch,
  data = adjusted_input
)

section("Example 3A output: dose effect adjusted for batch")
print(adjusted_fit$coefficient_table, digits = 12)

section("Example 3B output: effect of a two-dose-step increase")
two_step <- bb_contrast(adjusted_fit, c(dose = 2))
print(two_step, digits = 12)

stopifnot(
  abs(coef(dose_fit)[["dose"]] - 0.35) < 0.001,
  abs(coef(adjusted_fit)[["dose"]] - 0.3) < 0.001,
  abs(coef(adjusted_fit)[["batchB"]] - 0.2) < 0.001,
  abs(two_step$estimate - 2 * coef(adjusted_fit)[["dose"]]) < 1e-12
)
