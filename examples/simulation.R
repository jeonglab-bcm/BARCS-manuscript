#!/usr/bin/env Rscript

# Reproducible demonstration of beta-binomial regression for a continuous
# sample-level phenotype.  Run from the repository root:
#
#   Rscript examples/simulation.R

source(file.path("R", "bbreg.R"))
source(file.path("R", "method_palette.R"))

set.seed(20260723)

n_genes <- 200L
guides_per_gene <- 3L
n_guides <- n_genes * guides_per_gene
n_null_genes <- 160L
dose_grid <- seq(-1.5, 1.5, length.out = 9L)
sample_data <- expand.grid(
  dose = dose_grid,
  batch = factor(sprintf("batch_%d", 1:3)),
  KEEP.OUT.ATTRS = FALSE
)
# MAGeCK-MLE treats the first design row as the baseline. Put a dose-zero,
# reference-batch sample first while retaining the full continuous design.
baseline_row <- sample_data$dose == 0 &
  sample_data$batch == "batch_1"
sample_data <- rbind(
  sample_data[baseline_row, , drop = FALSE],
  sample_data[!baseline_row, , drop = FALSE]
)
n_samples <- nrow(sample_data)
library_size <- round(runif(n_samples, min = 40000, max = 120000))

guide <- sprintf("sg%04d", seq_len(n_guides))
gene_id <- sprintf("gene%04d", seq_len(n_genes))
gene <- rep(gene_id, each = guides_per_gene)
true_gene_effect <- c(
  rep(0, n_null_genes),
  rep(0.9, (n_genes - n_null_genes) / 2),
  rep(-0.9, (n_genes - n_null_genes) / 2)
)
true_effect <- rep(true_gene_effect, each = guides_per_gene)
truth_gene <- true_gene_effect != 0
baseline <- runif(n_guides, min = -7.4, max = -6.5)
batch_effect <- c(batch_1 = 0, batch_2 = 0.20, batch_3 = -0.15)
true_rho <- 0.0015
beta_precision <- 1 / true_rho - 1

counts <- matrix(
  0,
  nrow = n_guides,
  ncol = n_samples,
  dimnames = list(guide, sprintf("sample_%02d", seq_len(n_samples)))
)

for (g in seq_len(n_guides)) {
  eta <- baseline[g] +
    true_effect[g] * sample_data$dose +
    unname(batch_effect[as.character(sample_data$batch)])
  mu <- plogis(eta)
  latent_p <- rbeta(
    n_samples,
    shape1 = mu * beta_precision,
    shape2 = (1 - mu) * beta_precision
  )
  counts[g, ] <- rbinom(n_samples, size = library_size, prob = latent_p)
}

bb_result <- bb_screen(
  counts = counts,
  data = sample_data,
  formula = ~ dose + batch,
  term = "dose",
  totals = library_size,
  guide = guide,
  gene = gene
)

# Run the official MAGeCK-MLE executable on the identical continuous design.
mageck_executable <- file.path(".venv", "bin", "mageck")
mageck_dir <- file.path("results", "simulation_mageck")
dir.create(mageck_dir, recursive = TRUE, showWarnings = FALSE)
mageck_count_path <- file.path(mageck_dir, "counts.tsv")
mageck_design_path <- file.path(mageck_dir, "design.tsv")
mageck_prefix <- file.path(mageck_dir, "mageck")
mageck_gene_path <- paste0(mageck_prefix, ".gene_summary.txt")

write.table(
  data.frame(sgRNA = guide, Gene = gene, counts, check.names = FALSE),
  mageck_count_path,
  sep = "\t", quote = FALSE, row.names = FALSE
)
mageck_design <- data.frame(
  samples = colnames(counts),
  Initial_condition = 1,
  dose = sample_data$dose,
  batch_2 = as.integer(sample_data$batch == "batch_2"),
  batch_3 = as.integer(sample_data$batch == "batch_3"),
  check.names = FALSE
)
write.table(
  mageck_design, mageck_design_path,
  sep = "\t", quote = FALSE, row.names = FALSE
)
if (!file.exists(mageck_executable)) {
  stop("Official MAGeCK is required at `.venv/bin/mageck`.", call. = FALSE)
}
mageck_status <- system2(
  mageck_executable,
  c(
    "mle",
    "-k", mageck_count_path,
    "-d", mageck_design_path,
    "-n", mageck_prefix,
    "--norm-method", "median",
    "--permutation-round", "1",
    "--threads", "4"
  )
)
if (mageck_status != 0 || !file.exists(mageck_gene_path)) {
  stop("Official MAGeCK-MLE simulation analysis failed.", call. = FALSE)
}
mageck_raw <- read.delim(
  mageck_gene_path, check.names = FALSE, stringsAsFactors = FALSE
)
mageck_index <- match(gene_id, mageck_raw$Gene)
mageck_result <- data.frame(
  estimate = mageck_raw[["dose|beta"]][mageck_index],
  std_error = abs(
    mageck_raw[["dose|beta"]][mageck_index] /
      mageck_raw[["dose|z"]][mageck_index]
  ),
  p_value = mageck_raw[["dose|wald-p-value"]][mageck_index],
  fdr = mageck_raw[["dose|wald-fdr"]][mageck_index]
)

# Deliberately misspecified comparison: binomial regression assumes that all
# variation is sequencing variation.  Its anti-conservative behavior shows why
# beta-binomial overdispersion is needed.
fit_naive <- function(k) {
  fit <- suppressWarnings(glm(
    cbind(k, library_size - k) ~ dose + batch,
    family = binomial(),
    data = sample_data
  ))
  coefficient <- summary(fit)$coefficients["dose", ]
  c(estimate = coefficient[["Estimate"]],
    std_error = coefficient[["Std. Error"]],
    p_value = coefficient[["Pr(>|z|)"]])
}
naive_result <- t(apply(counts, 1L, fit_naive))
naive_fdr <- p.adjust(naive_result[, "p_value"], method = "BH")

aggregate_gene_result <- function(estimate, p_value) {
  index <- split(seq_along(gene), gene)
  gene_estimate <- vapply(index, function(i) {
    median(estimate[i], na.rm = TRUE)
  }, numeric(1L))
  gene_p <- vapply(index, function(i) {
    keep <- is.finite(estimate[i]) & is.finite(p_value[i])
    if (!any(keep)) return(NA_real_)
    signed_z <- sign(estimate[i][keep]) *
      qnorm(
        pmax(p_value[i][keep] / 2, .Machine$double.xmin),
        lower.tail = FALSE
      )
    2 * pnorm(-abs(sum(signed_z) / sqrt(length(signed_z))))
  }, numeric(1L))
  combined_z <- qnorm(
    pmax(gene_p / 2, .Machine$double.xmin), lower.tail = FALSE
  )
  gene_se <- abs(gene_estimate / combined_z)
  data.frame(
    estimate = gene_estimate,
    std_error = gene_se,
    p_value = gene_p,
    fdr = p.adjust(gene_p, method = "BH")
  )
}

naive_gene <- aggregate_gene_result(
  naive_result[, "estimate"], naive_result[, "p_value"]
)
bb_gene <- aggregate_gene_result(bb_result$estimate, bb_result$p_value)
bb_moderated_result <- bb_moderate_dispersion(bb_result, trend = TRUE)
bb_moderated_gene <- aggregate_gene_result(
  bb_moderated_result$estimate,
  bb_moderated_result$p_value
)

summarize_method <- function(method, estimate, standard_error, p_value,
                             adjusted_p_value) {
  selected <- adjusted_p_value < 0.05
  data.frame(
    method = method,
    null_type_I_at_0.05 = mean(
      p_value[!truth_gene] < 0.05, na.rm = TRUE
    ),
    power_at_0.05 = mean(p_value[truth_gene] < 0.05, na.rm = TRUE),
    discoveries_at_FDR_0.05 = sum(selected, na.rm = TRUE),
    empirical_FDR = if (any(selected, na.rm = TRUE)) {
      mean(!truth_gene[selected], na.rm = TRUE)
    } else {
      0
    },
    effect_RMSE = sqrt(
      mean((estimate - true_gene_effect)^2, na.rm = TRUE)
    ),
    coverage_95 = mean(
      abs(estimate - true_gene_effect) <=
        qt(0.975, df = n_samples - 4L) * standard_error,
      na.rm = TRUE
    )
  )
}

simulation_summary <- rbind(
  summarize_method(
    "naive binomial z",
    naive_gene$estimate,
    naive_gene$std_error,
    naive_gene$p_value,
    naive_gene$fdr
  ),
  summarize_method(
    "beta-binomial t",
    bb_gene$estimate,
    bb_gene$std_error,
    bb_gene$p_value,
    bb_gene$fdr
  ),
  summarize_method(
    "beta-binomial moderated t",
    bb_moderated_gene$estimate,
    bb_moderated_gene$std_error,
    bb_moderated_gene$p_value,
    bb_moderated_gene$fdr
  ),
  summarize_method(
    "official MAGeCK-MLE Wald",
    mageck_result$estimate,
    mageck_result$std_error,
    mageck_result$p_value,
    mageck_result$fdr
  )
)

write.csv(
  simulation_summary,
  file.path("results", "simulation_summary.csv"),
  row.names = FALSE
)

null_guide <- rep(!truth_gene, each = guides_per_gene)
simulation_diagnostics <- data.frame(
  statistic = c(
    "guide_null_p_below_0_05",
    "null_guides_at_rho_lower_boundary",
    "all_guides_at_rho_lower_boundary"
  ),
  value = c(
    mean(bb_result$p_value[null_guide] < 0.05, na.rm = TRUE),
    mean(bb_result$rho[null_guide] == 0, na.rm = TRUE),
    mean(bb_result$rho == 0, na.rm = TRUE)
  )
)
write.csv(
  simulation_diagnostics,
  file.path("results", "simulation_diagnostics.csv"),
  row.names = FALSE
)

example_index <- which(true_effect != 0 & bb_result$fdr < 0.05)[1L]
example_fit <- bbreg(
  count = counts[example_index, ],
  total = library_size,
  formula = ~ dose + batch,
  data = sample_data
)
example_gene_index <- match(gene[example_index], gene_id)
example_table <- data.frame(
  guide = guide[example_index],
  true_effect = true_effect[example_index],
  bb_estimate = example_fit$coefficient_table["dose", "estimate"],
  bb_std_error = example_fit$coefficient_table["dose", "std_error"],
  bb_t_value = example_fit$coefficient_table["dose", "t_value"],
  bb_df = example_fit$coefficient_table["dose", "df"],
  bb_p_value = example_fit$coefficient_table["dose", "p_value"],
  estimated_rho = example_fit$rho,
  mageck_estimate = mageck_result$estimate[example_gene_index],
  mageck_std_error = mageck_result$std_error[example_gene_index],
  mageck_p_value = mageck_result$p_value[example_gene_index],
  row.names = NULL,
  check.names = FALSE
)
write.csv(
  example_table,
  file.path("results", "worked_example.csv"),
  row.names = FALSE
)

pdf(file.path("figures", "simulation_diagnostics.pdf"),
    width = 9, height = 4.5, onefile = TRUE)
old_par <- par(mfrow = c(1, 2), mar = c(4.2, 4.5, 2.2, 1))

null_bb_p <- sort(
  bb_gene$p_value[!truth_gene & is.finite(bb_gene$p_value)]
)
null_mageck_p <- sort(
  mageck_result$p_value[
    !truth_gene & is.finite(mageck_result$p_value)
  ]
)
null_naive_p <- sort(
  naive_gene$p_value[
    !truth_gene & is.finite(naive_gene$p_value)
  ]
)
expected <- ppoints(length(null_bb_p))
plot(
  expected, null_naive_p,
  pch = 16, cex = 0.55,
  col = adjustcolor(barcs_method_colours[["Naive binomial"]], alpha.f = 0.65),
  xlab = "Expected null p-value",
  ylab = "Observed null p-value",
  main = "(A) Null calibration"
)
points(
  expected, null_bb_p, pch = 16, cex = 0.55,
  col = adjustcolor(barcs_method_colours[["BARCS"]], alpha.f = 0.65)
)
points(
  ppoints(length(null_mageck_p)), null_mageck_p,
  pch = 16, cex = 0.55,
  col = adjustcolor(barcs_method_colours[["MAGeCK"]], alpha.f = 0.65)
)
abline(0, 1, lty = 2)
legend(
  "topleft",
  legend = c(
    "naive binomial z", "BARCS t", "official MAGeCK-MLE Wald"
  ),
  pch = 16,
  col = unname(barcs_method_colours[c(
    "Naive binomial", "BARCS", "MAGeCK"
  )]),
  bty = "n",
  cex = 0.78
)

dose_prediction <- seq(min(sample_data$dose), max(sample_data$dose),
                       length.out = 200)
prediction_data <- data.frame(
  dose = dose_prediction,
  batch = factor("batch_1", levels = levels(sample_data$batch))
)
prediction_matrix <- model.matrix(~ dose + batch, prediction_data)
prediction <- plogis(drop(prediction_matrix %*% coef(example_fit)))
observed <- counts[example_index, ] / library_size
plot(
  sample_data$dose, observed * 1e6,
  pch = 16,
  col = c("#4D4D4D", "#969696", "#CCCCCC")[
    as.integer(sample_data$batch)
  ],
  xlab = "Standardized continuous phenotype",
  ylab = "Guide abundance (CPM)",
  main = sprintf("(B) %s: true slope = %.1f", guide[example_index],
                 true_effect[example_index])
)
lines(
  dose_prediction, prediction * 1e6, lwd = 2,
  col = barcs_method_colours[["BARCS"]]
)
legend(
  "topleft",
  legend = c(levels(sample_data$batch), "fitted batch 1"),
  pch = c(16, 16, 16, NA), lty = c(NA, NA, NA, 1),
  col = c(
    "#4D4D4D", "#969696", "#CCCCCC",
    barcs_method_colours[["BARCS"]]
  ),
  bty = "n"
)

par(old_par)
dev.off()

print(simulation_summary, row.names = FALSE, digits = 3)
cat("\nWorked guide:\n")
print(example_table, row.names = FALSE, digits = 4)
