#!/usr/bin/env Rscript

# Truth-based CRISPulator FACS benchmark for three BARCS gene statistics.
#
# All methods receive the identical guide-level beta-binomial fits:
#   1. BARCS-original: historical signed-z guide aggregation;
#   2. BARCS-partial: random-effects partial pooling of guide coefficients;
#   3. BARCS-EB: empirical-Bayes moderation of guide heterogeneity.
#
# The benchmark deliberately excludes unrelated count methods. Its purpose is
# to isolate the consequence of changing only the guide-to-gene statistic.

options(stringsAsFactors = FALSE)
source(file.path("R", "method_palette.R"))
analysis_protocol <- "barcs-three-methods-v1"

data_dir <- Sys.getenv(
  "CRISPULATOR_DATA_DIR",
  file.path("data", "derived", "crispulator_facs")
)
result_dir <- Sys.getenv(
  "CRISPULATOR_RESULT_DIR",
  file.path("results", "crispulator_facs", "three_methods")
)
figure_path <- Sys.getenv(
  "CRISPULATOR_FIGURE_PATH",
  file.path("figures", "crispulator_facs_benchmark.pdf")
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(figure_path), recursive = TRUE, showWarnings = FALSE)

required <- file.path(
  data_dir,
  c("counts.tsv", "sample_design.tsv", "guide_truth.tsv", "gene_truth.tsv")
)
if (!all(file.exists(required))) {
  stop(
    "Run `julia --project=julia julia/simulate_crispulator_facs.jl` first.",
    call. = FALSE
  )
}

count_table <- read.delim(required[1L], check.names = FALSE)
sample_data <- read.delim(required[2L], check.names = FALSE)
guide_truth <- read.delim(required[3L], check.names = FALSE)
gene_truth <- read.delim(required[4L], check.names = FALSE)
gene_truth$active <- tolower(as.character(gene_truth$active)) == "true"
counts <- as.matrix(count_table[, sample_data$sample, drop = FALSE])
storage.mode(counts) <- "double"
stopifnot(
  identical(count_table$guide, guide_truth$guide),
  identical(count_table$gene, guide_truth$gene),
  all(colSums(counts) > 0),
  all(table(sample_data$replicate) == length(unique(sample_data$sample_type)))
)
n_replicates <- length(unique(sample_data$replicate))
if (n_replicates < 1L) {
  stop("At least one screen replicate is required.")
}

cpp_library <- file.path("CB2", "src", paste0("CB2", .Platform$dynlib.ext))
if (file.exists(cpp_library) &&
    requireNamespace("Rcpp", quietly = TRUE) &&
    requireNamespace("RcppArmadillo", quietly = TRUE)) {
  suppressPackageStartupMessages(library(Rcpp))
  suppressPackageStartupMessages(library(RcppArmadillo))
  if (!"CB2" %in% names(getLoadedDLLs())) {
    dyn.load(cpp_library)
  }
  source(file.path("CB2", "R", "RcppExports.R"))
}
source(file.path("R", "bbreg.R"))

conditional_normal_mean <- function(lower, upper) {
  mass <- pnorm(upper) - pnorm(lower)
  (dnorm(lower) - dnorm(upper)) / mass
}
q25 <- qnorm(0.25)
sample_data$phenotype_z <- NA_real_
sample_data$phenotype_z[sample_data$sample_type == "low"] <-
  conditional_normal_mean(-Inf, q25)
sample_data$phenotype_z[sample_data$sample_type == "high"] <-
  conditional_normal_mean(-q25, Inf)
sample_data$phenotype_z[sample_data$sample_type == "bulk"] <- 0
sample_data$replicate <- factor(sample_data$replicate)

method_labels <- c(
  original = "BARCS-original",
  partial = "BARCS-partial",
  eb = "BARCS-EB"
)
method_colours <- c(
  `BARCS-original` = "#0072B2",
  `BARCS-partial` = "#009E73",
  `BARCS-EB` = "#D55E00"
)

fit_three_methods <- function(label, sample_index,
                              term = "phenotype_z") {
  design <- droplevels(sample_data[sample_index, , drop = FALSE])
  count_subset <- counts[, sample_index, drop = FALSE]
  if (term == "phenotype_z") {
    formula <- if (n_replicates >= 2L) {
      ~ phenotype_z + replicate
    } else {
      ~ phenotype_z
    }
  } else {
    design$bulk_indicator <- as.integer(design$sample_type == "bulk")
    formula <- ~ bulk_indicator + replicate
  }

  fit_start <- proc.time()
  guide_result <- bb_screen(
    counts = count_subset,
    totals = colSums(count_subset),
    data = design,
    formula = formula,
    term = term,
    guide = guide_truth$guide,
    gene = guide_truth$gene,
    min_total_count = 30,
    ncores = as.integer(Sys.getenv("BARCS_NCORES", "4"))
  )
  guide_seconds <- unname((proc.time() - fit_start)[["elapsed"]])
  negative_control <- guide_truth$class == "negcontrol"
  calibrated <- bb_calibrate_controls(
    guide_result,
    negative_control,
    alpha = 0.05
  )

  method_start <- proc.time()
  original <- bb_gene_original(calibrated, min_guides = 1L)
  original_seconds <- unname((proc.time() - method_start)[["elapsed"]])
  method_start <- proc.time()
  partial <- bb_gene_partial_pool(
    guide_result,
    control = negative_control,
    min_guides = 2L,
    alpha = 0.05,
    min_control_genes = 10L
  )
  partial_seconds <- unname((proc.time() - method_start)[["elapsed"]])
  method_start <- proc.time()
  eb <- bb_gene_eb_moderate(
    guide_result,
    control = negative_control,
    min_guides = 2L,
    prior_df = 4,
    alpha = 0.05,
    min_control_genes = 10L
  )
  eb_seconds <- unname((proc.time() - method_start)[["elapsed"]])

  write.csv(
    guide_result,
    gzfile(file.path(result_dir, paste0(label, "_guide_results.csv.gz"))),
    row.names = FALSE
  )
  write.csv(
    calibrated,
    gzfile(file.path(
      result_dir, paste0(label, "_calibrated_guide_results.csv.gz")
    )),
    row.names = FALSE
  )
  results <- list(original = original, partial = partial, eb = eb)
  for (method in names(results)) {
    write.csv(
      results[[method]],
      file.path(
        result_dir,
        paste0(label, "_", method, "_gene_results.csv")
      ),
      row.names = FALSE
    )
  }
  write.csv(
    data.frame(
      analysis_protocol = analysis_protocol,
      design = label,
      method = unname(method_labels),
      shared_guide_fit_seconds = guide_seconds,
      gene_aggregation_seconds = c(
        original_seconds, partial_seconds, eb_seconds
      ),
      converged_fraction = mean(guide_result$converged),
      guide_control_scale = attr(calibrated, "control_scale"),
      gene_null_scale = c(
        NA_real_,
        attr(partial, "null_scale"),
        attr(eb, "null_scale")
      ),
      prior_tau2 = c(NA_real_, NA_real_, attr(eb, "prior_tau2")),
      prior_df = c(NA_real_, NA_real_, attr(eb, "prior_df"))
    ),
    file.path(result_dir, paste0(label, "_runtime.csv")),
    row.names = FALSE
  )
  results
}

three_sample_index <- sample_data$sample_type %in% c("low", "bulk", "high")
tail_index <- sample_data$sample_type %in% c("low", "high")
bulk_index <- sample_data$sample_type %in% c("input", "bulk")

design_results <- list(
  `Low + bulk + high` = fit_three_methods(
    "low_bulk_high", three_sample_index
  )
)
if (n_replicates >= 2L) {
  design_results[["Two 25% tails"]] <- fit_three_methods(
    "two_tails", tail_index
  )
  design_results[["Unsorted 0-100%"]] <- fit_three_methods(
    "bulk_vs_input", bulk_index, term = "bulk_indicator"
  )
}

auroc <- function(truth, score) {
  score_rank <- rank(score, ties.method = "average")
  n_positive <- sum(truth)
  n_negative <- sum(!truth)
  (
    sum(score_rank[truth]) -
      n_positive * (n_positive + 1) / 2
  ) / (n_positive * n_negative)
}

average_precision <- function(truth, score) {
  ordered_truth <- truth[order(score, decreasing = TRUE)]
  precision_at_rank <- cumsum(ordered_truth) / seq_along(ordered_truth)
  mean(precision_at_rank[ordered_truth])
}

evaluate <- function(method, gene_result, design) {
  assessed <- merge(
    gene_truth,
    gene_result[, c("gene", "estimate", "p_value", "fdr")],
    by = "gene",
    all.x = TRUE
  )
  assessed <- assessed[
    is.finite(assessed$estimate) &
      is.finite(assessed$p_value) &
      is.finite(assessed$fdr),
    ,
    drop = FALSE
  ]
  score <- -log10(pmax(assessed$p_value, .Machine$double.xmin))
  called <- assessed$fdr < 0.10
  sign_match <- sign(assessed$estimate) == assessed$expected_sign
  active <- assessed$active
  true_positive <- sum(called & active)
  false_positive <- sum(called & !active)
  false_negative <- sum(!called & active)
  precision <- if ((true_positive + false_positive) > 0) {
    true_positive / (true_positive + false_positive)
  } else {
    0
  }
  recall <- if ((true_positive + false_negative) > 0) {
    true_positive / (true_positive + false_negative)
  } else {
    0
  }
  f1 <- if ((precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    0
  }
  negative_control <- assessed$class == "negcontrol"
  data.frame(
    analysis_protocol = analysis_protocol,
    method = method,
    design = design,
    genes = nrow(assessed),
    active_genes = sum(active),
    auroc = auroc(active, score),
    average_precision = average_precision(active, score),
    effect_spearman_active = cor(
      assessed$estimate[active],
      assessed$theoretical_phenotype[active],
      method = "spearman"
    ),
    direction_accuracy_active = mean(sign_match[active]),
    calls_fdr_0_10 = sum(called),
    directional_recall_fdr_0_10 =
      mean(called[active] & sign_match[active]),
    empirical_fdp_fdr_0_10 =
      if (sum(called) > 0) false_positive / sum(called) else 0,
    f1_fdr_0_10 = f1,
    negative_control_p_below_0_05 =
      mean(assessed$p_value[negative_control] < 0.05)
  )
}

metric_rows <- list()
effect_rows <- list()
row_index <- 0L
for (design_name in names(design_results)) {
  for (method_key in names(method_labels)) {
    row_index <- row_index + 1L
    gene_result <- design_results[[design_name]][[method_key]]
    metric_rows[[row_index]] <- evaluate(
      method_labels[[method_key]], gene_result, design_name
    )
    effect <- merge(
      gene_truth[
        , c("gene", "class", "active", "theoretical_phenotype")
      ],
      gene_result[, c("gene", "estimate", "p_value", "fdr")],
      by = "gene"
    )
    effect$method <- method_labels[[method_key]]
    effect$design <- design_name
    effect_rows[[row_index]] <- effect
  }
}
metrics <- do.call(rbind, metric_rows)
effects <- do.call(rbind, effect_rows)
write.csv(
  metrics,
  file.path(result_dir, "benchmark_metrics.csv"),
  row.names = FALSE
)
write.csv(
  effects,
  file.path(result_dir, "gene_effects_with_truth.csv"),
  row.names = FALSE
)

primary <- metrics[metrics$design == "Low + bulk + high", ]
pdf(figure_path, width = 10.5, height = 4.2, useDingbats = FALSE)
layout(matrix(1:3, nrow = 1))

par(mar = c(7, 4.2, 2.6, 0.8))
metric_matrix <- rbind(
  AUROC = primary$auroc,
  `Average precision` = primary$average_precision,
  `Directional recall` = primary$directional_recall_fdr_0_10
)
matplot(
  seq_len(nrow(metric_matrix)),
  metric_matrix,
  type = "b",
  pch = 16,
  lty = 1,
  lwd = 2,
  col = unname(method_colours[primary$method]),
  ylim = c(0, 1),
  xaxt = "n",
  xlab = "",
  ylab = "Truth-recovery metric",
  main = "A  Low + bulk + high",
  bty = "l"
)
axis(
  1,
  at = seq_len(nrow(metric_matrix)),
  labels = rownames(metric_matrix),
  las = 2,
  cex.axis = 0.72
)
legend(
  "bottomleft",
  legend = primary$method,
  col = unname(method_colours[primary$method]),
  lty = 1,
  lwd = 2,
  pch = 16,
  bty = "n",
  cex = 0.72
)

par(mar = c(7, 4.2, 2.6, 0.8))
fdp_positions <- barplot(
  primary$empirical_fdp_fdr_0_10,
  col = unname(method_colours[primary$method]),
  border = NA,
  ylim = c(0, max(0.12, primary$empirical_fdp_fdr_0_10)),
  names.arg = primary$method,
  las = 2,
  cex.names = 0.72,
  ylab = "Realized FDP",
  main = "B  FDR 0.10 calls"
)
abline(h = 0.10, lty = 2, col = "#555555")
text(
  fdp_positions,
  primary$empirical_fdp_fdr_0_10,
  labels = sprintf("%.3f", primary$empirical_fdp_fdr_0_10),
  pos = 3,
  cex = 0.72
)

par(mar = c(7, 4.2, 2.6, 0.8))
null_positions <- barplot(
  primary$negative_control_p_below_0_05,
  col = unname(method_colours[primary$method]),
  border = NA,
  ylim = c(0, max(0.08, primary$negative_control_p_below_0_05)),
  names.arg = primary$method,
  las = 2,
  cex.names = 0.72,
  ylab = "Negative-control fraction",
  main = "C  Null p < 0.05"
)
abline(h = 0.05, lty = 2, col = "#555555")
text(
  null_positions,
  pmax(primary$negative_control_p_below_0_05, 0.002),
  labels = sprintf(
    "%.3f", primary$negative_control_p_below_0_05
  ),
  pos = 3,
  cex = 0.72
)
dev.off()

print(metrics)
