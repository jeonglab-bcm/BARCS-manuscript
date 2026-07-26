#!/usr/bin/env Rscript

source(file.path("R", "bbreg.R"))

assert_close <- function(actual, expected, tolerance = 1e-6, message = NULL) {
  if (any(!is.finite(c(actual, expected))) ||
      max(abs(actual - expected)) > tolerance) {
    stop(
      if (is.null(message)) {
        sprintf("Expected %s, got %s.", expected, actual)
      } else {
        message
      },
      call. = FALSE
    )
  }
}

set.seed(404)
data_equal <- data.frame(
  dose = rep(seq(-1, 1, length.out = 8), 2),
  batch = factor(rep(c("a", "b"), each = 8))
)
total_equal <- rep(50000, nrow(data_equal))
mu_equal <- plogis(-6.2 + 0.7 * data_equal$dose +
                     0.15 * (data_equal$batch == "b"))
count_equal <- rbinom(length(mu_equal), total_equal, mu_equal)

fit <- bbreg(count_equal, total_equal, ~ dose + batch, data_equal)
reference <- glm(
  cbind(count_equal, total_equal - count_equal) ~ dose + batch,
  family = binomial(),
  data = data_equal
)

# With equal library sizes, beta-binomial weights differ from binomial weights
# by a common multiplier.  The coefficient estimates must therefore agree.
assert_close(
  coef(fit),
  coef(reference),
  tolerance = 1e-5,
  message = "Equal-library beta-binomial and binomial coefficients disagree."
)
stopifnot(
  inherits(fit, "bbreg"),
  fit$converged,
  fit$df.residual == nrow(data_equal) - length(coef(fit)),
  fit$rho >= 0,
  fit$rho < 1,
  all(fitted(fit) > 0 & fitted(fit) < 1)
)

# A one-replicate low--bulk--high screen has three observations and two
# coefficients. It is estimable, but its single residual degree of freedom
# makes it a diagnostic/ranking case rather than a confirmatory design.
single_screen_data <- data.frame(
  phenotype_z = c(-1.271, 0, 1.271)
)
single_screen_total <- rep(50000, 3L)
single_screen_count <- c(70, 100, 145)
single_screen_fit <- bbreg(
  single_screen_count,
  single_screen_total,
  ~ phenotype_z,
  single_screen_data
)
stopifnot(
  single_screen_fit$converged,
  single_screen_fit$df.residual == 1L,
  is.finite(coef(single_screen_fit)[["phenotype_z"]]),
  is.finite(
    single_screen_fit$coefficient_table["phenotype_z", "p_value"]
  )
)

dose_contrast <- bb_contrast(fit, c(dose = 1))
assert_close(dose_contrast$estimate, coef(fit)[["dose"]])
assert_close(
  dose_contrast$p_value,
  fit$coefficient_table["dose", "p_value"]
)

# A genuine beta-binomial sample should produce a positive dispersion estimate
# and a larger standard error than a binomial-only fit.
set.seed(405)
rho_true <- 0.003
precision <- 1 / rho_true - 1
latent <- rbeta(
  length(mu_equal),
  mu_equal * precision,
  (1 - mu_equal) * precision
)
count_overdispersed <- rbinom(length(mu_equal), total_equal, latent)
fit_overdispersed <- bbreg(
  count_overdispersed, total_equal, ~ dose + batch, data_equal
)
reference_overdispersed <- glm(
  cbind(count_overdispersed, total_equal - count_overdispersed) ~ dose + batch,
  family = binomial(),
  data = data_equal
)
stopifnot(
  fit_overdispersed$rho > 0,
  fit_overdispersed$coefficient_table["dose", "std_error"] >
    summary(reference_overdispersed)$coefficients["dose", "Std. Error"]
)

screen_counts <- rbind(
  guide_a = count_equal,
  guide_b = count_overdispersed
)
screen <- bb_screen(
  screen_counts,
  data_equal,
  ~ dose + batch,
  term = "dose",
  totals = total_equal,
  gene = c("gene_a", "gene_b")
)
stopifnot(
  nrow(screen) == 2L,
  all(c("gene", "guide", "estimate", "std_error", "t_value", "df",
        "p_value", "rho", "fdr") %in% names(screen)),
  all(screen$fdr >= screen$p_value, na.rm = TRUE)
)

control_t <- 1.8 * qt((seq_len(100) - 0.5) / 100, df = 8)
control_result <- data.frame(
  estimate = rep(1, 101),
  std_error = rep(1, 101),
  t_value = c(control_t, 5),
  df = rep(8, 101),
  p_value = 2 * pt(-abs(c(control_t, 5)), df = 8)
)
control_result$fdr <- p.adjust(control_result$p_value, method = "BH")
calibrated <- bb_calibrate_controls(
  control_result,
  c(rep(TRUE, 100), FALSE)
)
stopifnot(
  attr(calibrated, "control_scale") > 1,
  identical(calibrated$raw_t_value, control_result$t_value),
  all(calibrated$p_value >= control_result$p_value),
  abs(mean(calibrated$p_value[seq_len(100)] < 0.05) - 0.05) <= 0.02
)

# Guide-consistency inference estimates one inverse-variance-weighted gene
# effect when sample-level degrees of freedom are too small for useful guide
# p-values. It does not combine guide p-values by Fisher or Stouffer.
set.seed(406)
consistency_input <- data.frame(
  gene = rep(c("control_a", "control_b", "null", "signal"), each = 5),
  estimate = c(
    rnorm(5, 0, 0.10),
    rnorm(5, 0, 0.10),
    rnorm(5, 0, 0.10),
    rnorm(5, 0.65, 0.08)
  ),
  std_error = rep(0.10, 20),
  converged = TRUE
)
consistency_control <- consistency_input$gene %in%
  c("control_a", "control_b")
consistency <- bb_gene_consistency(
  consistency_input,
  control = consistency_control,
  min_control_genes = 2
)
stopifnot(
  nrow(consistency) == 4L,
  all(c("std_error", "raw_statistic", "statistic",
        "guide_direction_agreement",
        "control_gene", "p_value", "fdr") %in% names(consistency)),
  attr(consistency, "null_scale") >= 1,
  attr(consistency, "control_genes") == 2L,
  consistency$p_value[consistency$gene == "signal"] <
    consistency$p_value[consistency$gene == "null"],
  consistency$guide_direction_agreement[
    consistency$gene == "signal"
  ] == 1
)

calibrated_consistency_input <- consistency_input
calibrated_consistency_input$raw_std_error <-
  calibrated_consistency_input$std_error
calibrated_consistency_input$std_error <-
  calibrated_consistency_input$std_error * 3
calibrated_consistency <- bb_gene_consistency(
  calibrated_consistency_input,
  control = consistency_control,
  min_control_genes = 2
)
assert_close(
  calibrated_consistency$raw_statistic,
  consistency$raw_statistic,
  message = "Gene consistency did not recover retained raw standard errors."
)

# The historical method must reproduce the signed-z calculation exactly,
# while the two guide-aware methods retain and moderate heterogeneity.
pooling_input <- data.frame(
  gene = rep(
    c("control_a", "control_b", "consistent", "discordant"),
    each = 5
  ),
  estimate = c(
    -0.04, 0.02, 0.01, -0.03, 0.04,
    0.03, -0.02, 0.04, -0.01, -0.03,
    -1.00, -0.90, -1.10, -1.00, -1.05,
    -1.00, -0.90, 0.80, 1.00, 0.20
  ),
  std_error = rep(0.20, 20),
  converged = TRUE
)
pooling_input$p_value <- 2 * pnorm(
  -abs(pooling_input$estimate / pooling_input$std_error)
)
pooling_control <- pooling_input$gene %in% c("control_a", "control_b")

original_gene <- bb_gene_original(pooling_input)
consistent_index <- pooling_input$gene == "consistent"
expected_original_z <- sum(
  pooling_input$estimate[consistent_index] /
    pooling_input$std_error[consistent_index]
) / sqrt(sum(consistent_index))
assert_close(
  original_gene$statistic[original_gene$gene == "consistent"],
  expected_original_z,
  message = "Original BARCS signed-z aggregation changed."
)

normal_gene <- bb_gene_normal(pooling_input, min_guides = 3L)
normal_consistent <- normal_gene[
  normal_gene$gene == "consistent", , drop = FALSE
]
expected_normal_beta <- pooling_input$estimate[consistent_index]
expected_normal_mean <- mean(expected_normal_beta)
expected_normal_sigma <- sd(expected_normal_beta)
expected_normal_standard_error <-
  expected_normal_sigma / sqrt(length(expected_normal_beta))
expected_normal_statistic <-
  expected_normal_mean / expected_normal_standard_error
assert_close(normal_consistent$estimate, expected_normal_mean)
assert_close(normal_consistent$sigma, expected_normal_sigma)
assert_close(
  normal_consistent$std_error,
  expected_normal_standard_error
)
assert_close(
  normal_consistent$statistic,
  expected_normal_statistic
)
assert_close(
  normal_consistent$student_p_value,
  2 * pt(
    -abs(expected_normal_statistic),
    df = length(expected_normal_beta) - 1L
  )
)
assert_close(
  normal_consistent$normal_p_value,
  2 * pnorm(-abs(expected_normal_statistic))
)
stopifnot(
  normal_consistent$df == length(expected_normal_beta) - 1L,
  normal_consistent$p_value == normal_consistent$student_p_value,
  normal_consistent$normal_p_value <
    normal_consistent$student_p_value,
  attr(normal_gene, "reference") == "student_t"
)
normal_gene_plugin <- bb_gene_normal(
  pooling_input,
  min_guides = 3L,
  reference = "normal"
)
normal_consistent_plugin <- normal_gene_plugin[
  normal_gene_plugin$gene == "consistent", , drop = FALSE
]
stopifnot(
  normal_consistent_plugin$p_value ==
    normal_consistent_plugin$normal_p_value,
  attr(normal_gene_plugin, "reference") == "normal"
)

partial_gene <- bb_gene_partial_pool(
  pooling_input,
  control = pooling_control,
  min_control_genes = 2
)
moderated_gene <- bb_gene_eb_moderate(
  pooling_input,
  control = pooling_control,
  min_control_genes = 2,
  prior_df = 4
)
partial_consistent <- partial_gene[
  partial_gene$gene == "consistent", , drop = FALSE
]
partial_discordant <- partial_gene[
  partial_gene$gene == "discordant", , drop = FALSE
]
stopifnot(
  partial_consistent$raw_tau2 < partial_discordant$raw_tau2,
  partial_consistent$guide_direction_agreement == 1,
  partial_discordant$guide_direction_agreement < 1,
  partial_consistent$p_value < partial_discordant$p_value,
  all(c(
    "tau2", "raw_tau2", "i_squared", "max_weight_fraction",
    "leave_one_out_max_change", "leave_one_out_sign_stable"
  ) %in% names(partial_gene)),
  attr(moderated_gene, "prior_df") == 4,
  is.finite(attr(moderated_gene, "prior_tau2")),
  all(moderated_gene$tau2 >= pmin(
    moderated_gene$raw_tau2,
    attr(moderated_gene, "prior_tau2")
  )),
  all(moderated_gene$tau2 <= pmax(
    moderated_gene$raw_tau2,
    attr(moderated_gene, "prior_tau2")
  ))
)

bad_response_failed <- inherits(
  try(bbreg(c(2, 3), c(1, 4), ~ 1, data.frame(x = 1:2)), silent = TRUE),
  "try-error"
)
bad_design_failed <- inherits(
  try(
    bbreg(
      count_equal,
      total_equal,
      ~ dose + duplicate,
      transform(data_equal, duplicate = dose)
    ),
    silent = TRUE
  ),
  "try-error"
)
stopifnot(bad_response_failed, bad_design_failed)

cat("All beta-binomial regression tests passed.\n")
