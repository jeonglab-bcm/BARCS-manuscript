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

section("Example 1: a two-group calculation that can be checked by hand")

p <- c(A = 0.0014, B = 0.0023)
v <- c(A = 4e-8, B = 5e-8)
x <- cbind(intercept = 1, group_B = c(0, 1))
precision <- 1 / v

information <- crossprod(x, precision * x)
score_target <- crossprod(x, precision * p)
beta <- drop(solve(information, score_target))
covariance <- solve(information)

direct <- c(
  effect = p[["B"]] - p[["A"]],
  std_error = sqrt(sum(v)),
  t_value = (p[["B"]] - p[["A"]]) / sqrt(sum(v))
)
regression <- c(
  effect = beta[[2]],
  std_error = sqrt(covariance[2, 2]),
  t_value = beta[[2]] / sqrt(covariance[2, 2])
)

print(data.frame(
  quantity = names(direct),
  direct_CB2 = unname(direct),
  regression = unname(regression)
), digits = 10, row.names = FALSE)
cat("Regression coefficients:\n")
print(setNames(beta, colnames(x)), digits = 10)
stopifnot(max(abs(direct - regression)) < 1e-12)

section("Example 2: stress-test the identity on 1,000 random problems")

set.seed(20260723)
random_error <- replicate(1000, {
  random_p <- runif(2, 1e-5, 0.05)
  variance_scale <- 10^runif(1, -10, -6)
  random_v <- variance_scale * runif(2, 0.5, 2)
  random_precision <- 1 / random_v
  random_information <- crossprod(x, random_precision * x)
  random_beta <- drop(solve(
    random_information,
    crossprod(x, random_precision * random_p)
  ))
  random_covariance <- solve(random_information)

  direct_effect <- random_p[2] - random_p[1]
  direct_se <- sqrt(sum(random_v))
  max(abs(c(
    random_beta[1] - random_p[1],
    random_beta[2] - direct_effect,
    random_covariance[2, 2] - sum(random_v),
    random_beta[2] / sqrt(random_covariance[2, 2]) -
      direct_effect / direct_se
  )))
})

cat(sprintf(
  "Largest discrepancy across 1,000 random problems: %.3g\n",
  max(random_error)
))
stopifnot(max(random_error) < 1e-8)

section("Example 3: use the original CB2 fit_ab() on read counts")

count_a <- matrix(c(74, 112, 91, 139), nrow = 1)
total_a <- matrix(c(52000, 81000, 69000, 97000), nrow = 1)
count_b <- matrix(c(128, 177, 164, 211), nrow = 1)
total_b <- matrix(c(55000, 76000, 72000, 93000), nrow = 1)

fit_a <- fit_ab(count_a, total_a)
fit_b <- fit_ab(count_b, total_b)
phat <- c(A = fit_a$phat, B = fit_b$phat)
vhat <- c(A = fit_a$vhat, B = fit_b$vhat)

legacy_information <- crossprod(x, (1 / vhat) * x)
legacy_beta <- drop(solve(
  legacy_information,
  crossprod(x, (1 / vhat) * phat)
))
legacy_covariance <- solve(legacy_information)

legacy_t <- diff(phat) / sqrt(sum(vhat))
regression_t <- legacy_beta[2] / sqrt(legacy_covariance[2, 2])
legacy_df <- sum(vhat)^2 /
  (vhat[1]^2 / 3 + vhat[2]^2 / 3)
legacy_p <- 2 * pt(-abs(legacy_t), df = legacy_df)
regression_p <- 2 * pt(-abs(regression_t), df = legacy_df)

comparison <- data.frame(
  quantity = c("effect", "standard error", "t", "df", "two-sided p"),
  original_CB2 = c(
    diff(phat),
    sqrt(sum(vhat)),
    legacy_t,
    legacy_df,
    legacy_p
  ),
  weighted_regression = c(
    legacy_beta[2],
    sqrt(legacy_covariance[2, 2]),
    regression_t,
    legacy_df,
    regression_p
  )
)
print(comparison, digits = 12, row.names = FALSE)
cat(sprintf(
  "Maximum absolute discrepancy: %.3g\n",
  max(abs(comparison$original_CB2 - comparison$weighted_regression))
))
stopifnot(
  max(abs(comparison$original_CB2 -
            comparison$weighted_regression)) < 1e-12
)

section("Example 4: show why default bbreg(~ group) is close, not identical")

sample_data <- data.frame(group = factor(rep(c("A", "B"), each = 4)))
all_count <- c(drop(count_a), drop(count_b))
all_total <- c(drop(total_a), drop(total_b))
logit_fit <- bbreg(
  count = all_count,
  total = all_total,
  formula = ~ group,
  data = sample_data
)

finite_sample <- data.frame(
  method = c("Original CB2", "Default CB2-Reg"),
  effect = c(diff(phat), coef(logit_fit)[["groupB"]]),
  effect_scale = c("raw proportion difference", "log-odds difference"),
  t_value = c(legacy_t, logit_fit$coefficient_table["groupB", "t_value"]),
  df = c(legacy_df, logit_fit$coefficient_table["groupB", "df"]),
  p_value = c(legacy_p, logit_fit$coefficient_table["groupB", "p_value"])
)
print(finite_sample, digits = 10, row.names = FALSE)
cat(
  "\nThe effects use different units, so compare the t-statistics rather",
  "than the raw effect numbers.\n"
)

section("Example 5: move the two groups closer and watch the tests converge")

p0 <- 0.002
delta <- c(8e-4, 4e-4, 2e-4, 1e-4, 5e-5, 2.5e-5)
var_a <- 1.2e-8
var_b <- 1.7e-8
local_equivalence <- do.call(rbind, lapply(delta, function(d) {
  p_a <- p0 - d / 2
  p_b <- p0 + d / 2
  t_raw <- d / sqrt(var_a + var_b)
  derivative_a <- 1 / (p_a * (1 - p_a))
  derivative_b <- 1 / (p_b * (1 - p_b))
  t_logit <- (qlogis(p_b) - qlogis(p_a)) /
    sqrt(derivative_a^2 * var_a + derivative_b^2 * var_b)
  data.frame(
    group_difference = d,
    raw_t = t_raw,
    logit_t = t_logit,
    logit_over_raw = t_logit / t_raw,
    absolute_gap = abs(t_logit - t_raw)
  )
}))
print(local_equivalence, digits = 8, row.names = FALSE)
stopifnot(
  tail(local_equivalence$absolute_gap, 1) <
    local_equivalence$absolute_gap[1]
)

section("Example 6: replace the group column with a continuous dose column")

dose <- rep(0:4, each = 2)
dose_total <- rep(c(80000, 100000), 5)
dose_probability <- plogis(-7 + 0.35 * dose)
dose_count <- round(dose_total * dose_probability)
dose_data <- data.frame(dose = dose)

cat("The design matrix is now:\n")
print(model.matrix(~ dose, dose_data))

dose_fit <- bbreg(
  count = dose_count,
  total = dose_total,
  formula = ~ dose,
  data = dose_data
)

cat("\nObserved guide proportions:\n")
print(data.frame(
  dose = dose,
  total = dose_total,
  count = dose_count,
  observed_proportion = dose_count / dose_total
), digits = 8, row.names = FALSE)

cat("\nFitted coefficients:\n")
print(dose_fit$coefficient_table, digits = 10)
cat(
  "\nThe fitted dose coefficient is",
  sprintf("%.6f", coef(dose_fit)[["dose"]]),
  "log-odds units per dose step; the data were created with 0.35.\n"
)

section("Conclusion")
cat(
  "Examples 1-3 show exact equality for the legacy two-group statistic.\n",
  "Example 4 shows the finite-sample difference made by the default logit ",
  "model.\n",
  "Example 5 shows the local convergence.\n",
  "Example 6 shows the generalization: the same second design-matrix column ",
  "can contain dose rather than only 0/1 group labels.\n",
  sep = ""
)
