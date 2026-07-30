#!/usr/bin/env Rscript

# Null beta-binomial calibration grid for the Supplement. The main grid varies
# the number of independent libraries, guide abundance, overdispersion, and
# guides per gene while keeping every tested group coefficient exactly zero.
# A focused Gaussian-copula arm adds within-gene guide dependence. A
# deterministic split of the null genes also tests whether calibration at the
# gene-aggregation level repairs error introduced by the Stouffer combiner.

options(stringsAsFactors = FALSE)

if (requireNamespace("BARCS", quietly = TRUE)) {
  suppressPackageStartupMessages(library(BARCS))
}
source(file.path("R", "bbreg.R"))

output_dir <- file.path("data", "derived")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sample_sizes <- c(4L, 6L, 8L, 12L)
mean_counts <- c(10, 100, 1000)
overdispersion_factors <- c(0, 5, 60, 180)
guides_per_gene_values <- c(3L, 5L)
seeds <- 3101:3105
library_total <- 100000L
n_genes <- 100L

simulate_counts <- function(n_samples, mean_count, rho, seed,
                            guides_per_gene, guide_correlation) {
  set.seed(seed)
  n_guides <- n_genes * guides_per_gene
  mu <- mean_count / library_total
  probability <- if (rho == 0) {
    matrix(mu, nrow = n_guides, ncol = n_samples)
  } else if (guide_correlation == 0) {
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
  } else {
    kappa <- 1 / rho - 1
    probability <- matrix(NA_real_, nrow = n_guides, ncol = n_samples)
    for (sample_index in seq_len(n_samples)) {
      shared_gene <- rep(rnorm(n_genes), each = guides_per_gene)
      guide_noise <- rnorm(n_guides)
      latent_normal <-
        sqrt(guide_correlation) * shared_gene +
        sqrt(1 - guide_correlation) * guide_noise
      uniform_score <- pmin(
        pmax(pnorm(latent_normal), .Machine$double.eps),
        1 - .Machine$double.eps
      )
      probability[, sample_index] <- qbeta(
        uniform_score,
        shape1 = mu * kappa,
        shape2 = (1 - mu) * kappa
      )
    }
    probability
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

split_calibrate_genes <- function(gene_result, alpha = 0.05) {
  gene_number <- as.integer(sub("^G", "", gene_result$gene))
  calibration <- gene_number %% 2L == 1L
  evaluation <- !calibration
  signed_z <- sign(gene_result$estimate) * qnorm(
    pmax(gene_result$p_value / 2, .Machine$double.xmin),
    lower.tail = FALSE
  )
  scale <- max(
    1,
    unname(quantile(
      abs(signed_z[calibration & is.finite(signed_z)]),
      1 - alpha,
      type = 8
    )) / qnorm(1 - alpha / 2)
  )
  calibrated_p <- 2 * pnorm(-abs(signed_z / scale))
  list(
    scale = scale,
    heldout_tests = sum(evaluation & is.finite(calibrated_p)),
    heldout_type_i = mean(
      calibrated_p[evaluation & is.finite(calibrated_p)] < alpha
    )
  )
}

summarize_fit <- function(result, moderated, scenario) {
  gene_result <- bb_gene_original(result)
  gene_calibration <- split_calibrate_genes(gene_result)
  data.frame(
    scenario,
    moderation = moderated,
    guide_tests = sum(is.finite(result$p_value)),
    guide_type_i = mean(result$p_value < 0.05, na.rm = TRUE),
    guide_any_bh_call = as.integer(any(result$fdr < 0.05, na.rm = TRUE)),
    gene_tests = sum(is.finite(gene_result$p_value)),
    gene_type_i = mean(gene_result$p_value < 0.05, na.rm = TRUE),
    gene_any_bh_call = as.integer(any(gene_result$fdr < 0.05, na.rm = TRUE)),
    gene_split_control_scale = gene_calibration$scale,
    gene_split_heldout_tests = gene_calibration$heldout_tests,
    gene_split_heldout_type_i = gene_calibration$heldout_type_i,
    rho_zero_fraction = mean(result$rho == 0, na.rm = TRUE),
    convergence_fraction = mean(result$converged, na.rm = TRUE)
  )
}

independent_grid <- expand.grid(
  samples = sample_sizes,
  mean_count = mean_counts,
  overdispersion_factor = overdispersion_factors,
  guides_per_gene = guides_per_gene_values,
  guide_correlation = 0,
  seed = seeds,
  KEEP.OUT.ATTRS = FALSE
)
# At mean count 10, the two highest dispersion factors produce mostly
# all-zero guide trajectories and too few finite fits for the moderation
# estimator. They are non-identifiable rather than calibration experiments.
independent_grid <- independent_grid[
  !(independent_grid$mean_count == 10 &
      independent_grid$overdispersion_factor >= 60),
  ,
  drop = FALSE
]
correlated_grid <- expand.grid(
  samples = 8L,
  mean_count = 100,
  overdispersion_factor = c(60, 180),
  guides_per_gene = guides_per_gene_values,
  guide_correlation = 0.4,
  seed = seeds,
  KEEP.OUT.ATTRS = FALSE
)
grid <- rbind(independent_grid, correlated_grid)
grid$rho <- grid$overdispersion_factor / (library_total - 1)

pieces <- vector("list", nrow(grid) * 2L)
piece_index <- 0L
for (row_index in seq_len(nrow(grid))) {
  scenario <- grid[row_index, ]
  message(
    "m=", scenario$samples,
    ", mean=", scenario$mean_count,
    ", (N-1)rho=", scenario$overdispersion_factor,
    ", guides=", scenario$guides_per_gene,
    ", guide-correlation=", scenario$guide_correlation,
    ", seed=", scenario$seed
  )
  group <- rep(c(0, 1), each = scenario$samples / 2)
  counts <- simulate_counts(
    scenario$samples,
    scenario$mean_count,
    scenario$rho,
    scenario$seed,
    scenario$guides_per_gene,
    scenario$guide_correlation
  )
  n_guides <- nrow(counts)
  guide <- sprintf("g%04d", seq_len(n_guides))
  gene <- rep(
    sprintf("G%03d", seq_len(n_genes)),
    each = scenario$guides_per_gene
  )
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
    library_total = library_total,
    overdispersion_factor = scenario$overdispersion_factor,
    rho = scenario$rho,
    guides_per_gene = scenario$guides_per_gene,
    guide_correlation = scenario$guide_correlation,
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
  "gene_any_bh_call", "gene_split_control_scale",
  "gene_split_heldout_type_i", "rho_zero_fraction",
  "convergence_fraction"
)
group_names <- c(
  "samples", "mean_count", "library_total", "overdispersion_factor",
  "rho", "guides_per_gene", "guide_correlation", "moderation"
)
summary_mean <- aggregate(
  replicate_results[, metric_names],
  by = replicate_results[, group_names],
  FUN = mean,
  na.rm = TRUE
)
summary_sd <- aggregate(
  replicate_results[, metric_names],
  by = replicate_results[, group_names],
  FUN = sd,
  na.rm = TRUE
)
names(summary_sd)[-(seq_along(group_names))] <- paste0(
  names(summary_sd)[-(seq_along(group_names))], "_sd"
)
summary_result <- merge(
  summary_mean,
  summary_sd,
  by = group_names,
  sort = TRUE
)
write.csv(
  summary_result,
  file.path(output_dir, "barcs_null_calibration_grid_summary.csv"),
  row.names = FALSE
)

print(summary_result, row.names = FALSE)
