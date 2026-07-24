#!/usr/bin/env Rscript

# Executable verification of two connections between original CB2 and CB2-Reg.
# The first is an exact legacy summary-scale representation, not strict nesting
# of the default estimator. The second is a local/asymptotic statement.

if (file.exists(file.path("CB2", "DESCRIPTION")) &&
    requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all("CB2", quiet = TRUE)
} else if (!requireNamespace("CB2", quietly = TRUE)) {
  stop("Install CB2, or run this script from the CB2-Reg repository root.")
}

count_a <- matrix(c(74, 112, 91, 139), nrow = 1)
total_a <- matrix(c(52000, 81000, 69000, 97000), nrow = 1)
count_b <- matrix(c(128, 177, 164, 211), nrow = 1)
total_b <- matrix(c(55000, 76000, 72000, 93000), nrow = 1)

fit_a <- fit_ab(count_a, total_a)
fit_b <- fit_ab(count_b, total_b)
phat <- c(A = fit_a$phat, B = fit_b$phat)
vhat <- c(A = fit_a$vhat, B = fit_b$vhat)

if (any(!is.finite(vhat)) || any(vhat <= 0)) {
  stop("The deterministic example must have two positive group variances.")
}

# The two rows are the sufficient group summaries. This is the compressed form
# of the sample-level GLS construction proved in the manuscript.
x <- cbind(intercept = 1, group_B = c(0, 1))
omega <- 1 / vhat
information <- crossprod(x, omega * x)
score_target <- crossprod(x, omega * phat)
beta_gls <- drop(solve(information, score_target))
cov_gls <- solve(information)

t_cb2 <- (phat[["B"]] - phat[["A"]]) / sqrt(sum(vhat))
t_gls <- beta_gls[[2]] / sqrt(cov_gls[2, 2])
n_a <- ncol(count_a)
n_b <- ncol(count_b)
df_cb2 <- sum(vhat)^2 /
  (vhat[["A"]]^2 / (n_a - 1) + vhat[["B"]]^2 / (n_b - 1))
p_cb2 <- 2 * pt(-abs(t_cb2), df = df_cb2)
p_gls <- 2 * pt(-abs(t_gls), df = df_cb2)

exact_checks <- c(
  intercept = beta_gls[[1]] - phat[["A"]],
  slope = beta_gls[[2]] - (phat[["B"]] - phat[["A"]]),
  slope_variance = cov_gls[2, 2] - sum(vhat),
  t_statistic = t_gls - t_cb2,
  p_value = p_gls - p_cb2
)
stopifnot(max(abs(exact_checks)) < 1e-12)

cat("Legacy summary-scale GLS representation\n")
print(data.frame(
  quantity = c("effect", "standard error", "t", "df", "two-sided p"),
  original_CB2 = c(
    phat[["B"]] - phat[["A"]],
    sqrt(sum(vhat)),
    t_cb2,
    df_cb2,
    p_cb2
  ),
  GLS_contrast = c(
    beta_gls[[2]],
    sqrt(cov_gls[2, 2]),
    t_gls,
    df_cb2,
    p_gls
  ),
  check.names = FALSE
), digits = 12, row.names = FALSE)
cat(sprintf("Maximum absolute discrepancy: %.3g\n\n", max(abs(exact_checks))))

# As the two probabilities approach a shared interior value, the delta-method
# logit statistic approaches the raw-proportion statistic. The variances are
# fixed here only to isolate the link transformation; the proposition supplies
# the sampling-asymptotic statement.
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
    difference = d,
    raw_t = t_raw,
    logit_t = t_logit,
    ratio = t_logit / t_raw,
    absolute_gap = abs(t_logit - t_raw)
  )
}))

cat("Local equivalence of raw- and logit-scale statistics\n")
print(local_equivalence, digits = 8, row.names = FALSE)
stopifnot(tail(local_equivalence$absolute_gap, 1) <
            local_equivalence$absolute_gap[1])
