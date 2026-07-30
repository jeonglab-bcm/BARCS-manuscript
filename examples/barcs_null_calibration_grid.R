#!/usr/bin/env Rscript

# Null beta-binomial calibration grid for the Supplement.  The grid varies the
# number of independent libraries, guide abundance, and intraclass correlation
# while keeping every tested group coefficient exactly zero.

options(stringsAsFactors = FALSE)

if (requireNamespace("BARCS", quietly = TRUE)) {
  suppressPackageStartupMessages(library(BARCS))
}
source(file.path("R", "bbreg.R"))

output_dir <- file.path("data", "derived")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sample_sizes <- c(4L, 6L, 8L, 12L)
mean_counts <- c(10, 100, 1000)
rho_values <- c(0, 5e-7, 5e-6)
seeds <- 3101:3105
library_total <- 1000000L
n_genes <- 100L
guides_per_gene <- 5L
n_guides <- n_genes * guides_per_gene

simulate_counts <- function(n_samples, mean_count, rho, seed) {
  set.seed(seed)
  mu <- mean_count / library_total
  probability <- if (rho == 0) {
    matrix(mu, nrow = n_guides, ncol = n_samples)
  } else {
    kappa <- 1 / rho - 1
    matrix(
      rbeta(
        n_guides * n_samples,
        shape1 = mu * kappa,
        shape2 = (1 - mu) * kappa
      ),
      nrow = n_guides,
      ncol = n_samples
    )
  }
  matrix(
    rbinom(
      n_guides * n_samples,
      size = library_total,
      prob = as.vector(probability)
    ),
    nrow = n_guides,
    ncol = n_samples
  )
}

summarize_fit <- function(result, moderated, scenario) {
  gene_result <- bb_gene_original(result)
  data.frame(
    scenario,
    moderation = moderated,
    guide_tests = sum(is.finite(result$p_value)),
    guide_type_i = mean(result$p_value < 0.05, na.rm = TRUE),
    guide_any_bh_call = as.integer(any(result$fdr < 0.05, na.rm = TRUE)),
    gene_tests = sum(is.finite(gene_result$p_value)),
    gene_type_i = mean(gene_result$p_value < 0.05, na.rm = TRUE),
    gene_any_bh_call = as.integer(any(gene_result$fdr < 0.05, na.rm = TRUE)),
    rho_zero_fraction = mean(result$rho == 0, na.rm = TRUE),
    convergence_fraction = mean(result$converged, na.rm = TRUE)
  )
}

grid <- expand.grid(
  samples = sample_sizes,
  mean_count = mean_counts,
  rho = rho_values,
  seed = seeds,
  KEEP.OUT.ATTRS = FALSE
)

pieces <- vector("list", nrow(grid) * 2L)
piece_index <- 0L
for (row_index in seq_len(nrow(grid))) {
  scenario <- grid[row_index, ]
  message(
    "m=", scenario$samples,
    ", mean=", scenario$mean_count,
    ", rho=", format(scenario$rho, scientific = TRUE),
    ", seed=", scenario$seed
  )
  group <- rep(c(0, 1), each = scenario$samples / 2)
  counts <- simulate_counts(
    scenario$samples,
    scenario$mean_count,
    scenario$rho,
    scenario$seed
  )
  guide <- sprintf("g%04d", seq_len(n_guides))
  gene <- rep(sprintf("G%03d", seq_len(n_genes)), each = guides_per_gene)
  fit <- bb_screen(
    counts = counts,
    totals = rep(library_total, scenario$samples),
    data = data.frame(group = group),
    formula = ~ group,
    term = "group",
    guide = guide,
    gene = gene,
    min_total_count = 1,
    ncores = min(4L, parallel::detectCores(logical = FALSE))
  )
  moderated <- bb_moderate_dispersion(
    fit,
    trend = TRUE,
    one_way = FALSE,
    borrow_df = TRUE
  )
  scenario_fields <- data.frame(
    samples = scenario$samples,
    mean_count = scenario$mean_count,
    rho = scenario$rho,
    seed = scenario$seed
  )
  piece_index <- piece_index + 1L
  pieces[[piece_index]] <- summarize_fit(fit, FALSE, scenario_fields)
  piece_index <- piece_index + 1L
  pieces[[piece_index]] <- summarize_fit(
    moderated, TRUE, scenario_fields
  )
}

replicate_results <- do.call(rbind, pieces)
write.csv(
  replicate_results,
  file.path(output_dir, "barcs_null_calibration_grid_replicates.csv"),
  row.names = FALSE
)

metric_names <- c(
  "guide_type_i", "guide_any_bh_call", "gene_type_i",
  "gene_any_bh_call", "rho_zero_fraction", "convergence_fraction"
)
summary_mean <- aggregate(
  replicate_results[, metric_names],
  by = replicate_results[, c(
    "samples", "mean_count", "rho", "moderation"
  )],
  FUN = mean,
  na.rm = TRUE
)
summary_sd <- aggregate(
  replicate_results[, metric_names],
  by = replicate_results[, c(
    "samples", "mean_count", "rho", "moderation"
  )],
  FUN = sd,
  na.rm = TRUE
)
names(summary_sd)[-(1:4)] <- paste0(
  names(summary_sd)[-(1:4)], "_sd"
)
summary_result <- merge(
  summary_mean,
  summary_sd,
  by = c("samples", "mean_count", "rho", "moderation"),
  sort = TRUE
)
write.csv(
  summary_result,
  file.path(output_dir, "barcs_null_calibration_grid_summary.csv"),
  row.names = FALSE
)

print(summary_result, row.names = FALSE)
