#!/usr/bin/env Rscript

# Waterbear IL2RA FACS data benchmark for four BARCS gene statistics.
#
# The same 12 low-coverage/high-MOI count columns, sample design, guide fits,
# and negative-control definition are used for:
#   1. BARCS-original;
#   2. BARCS-NORM;
#   3. BARCS-partial;
#   4. BARCS-EB.
#
# The 26 directionally validated genes and the selected 33-gene follow-up
# panel are supporting truth sets. The latter is not an unbiased genome-wide
# negative set because candidates were selected for experimental follow-up.

options(stringsAsFactors = FALSE)
source(file.path("R", "method_palette.R"))
analysis_protocol <- "barcs-four-methods-v1"

raw_dir <- file.path("data", "raw", "waterbear")
result_dir <- file.path("results", "waterbear_facs", "three_methods")
figure_path <- file.path("figures", "waterbear_facs_benchmark.pdf")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(figure_path), recursive = TRUE, showWarnings = FALSE)

count_source <- file.path(
  raw_dir,
  "GSE242880_high_low_MOI_screens_D1_D2_D3_2022-04-23.count.txt.gz"
)
validation_source <- file.path(
  raw_dir, "cd25_validated_targets_revised.csv"
)
if (!all(file.exists(c(count_source, validation_source)))) {
  stop(
    "Waterbear processed counts and validation table are required under ",
    "`data/raw/waterbear`.",
    call. = FALSE
  )
}

# Use BARCS's own compiled RcppArmadillo kernels when the package is installed.
# The base-R fallback in R/bbreg.R remains valid but is slower.
if (requireNamespace("BARCS", quietly = TRUE)) {
  suppressPackageStartupMessages(library(BARCS))
}
source(file.path("R", "bbreg.R"))

raw <- read.delim(gzfile(count_source), check.names = FALSE)
sample_columns <- grep(
  "^Low_Coverage_High_MOI_D[123]_Q[1-4]$",
  names(raw),
  value = TRUE
)
if (length(sample_columns) != 12L) {
  stop("Expected four low-coverage/high-MOI bins for each of three donors.")
}
counts <- as.matrix(raw[, sample_columns, drop = FALSE])
storage.mode(counts) <- "double"
sample_data <- data.frame(
  sample = sample_columns,
  donor = factor(sub(".*_(D[123])_Q[1-4]$", "\\1", sample_columns)),
  bin = as.integer(sub(".*_Q([1-4])$", "\\1", sample_columns))
)

bin_mass <- c(0.2, 0.3, 0.3, 0.2)
bin_cut <- qnorm(c(0, cumsum(bin_mass)))
bin_location <- vapply(seq_along(bin_mass), function(index) {
  (dnorm(bin_cut[index]) - dnorm(bin_cut[index + 1L])) /
    bin_mass[index]
}, numeric(1L))
sample_data$marker_z <- bin_location[sample_data$bin]

write.table(
  data.frame(
    sgRNA = raw$sgRNA,
    Gene = raw$Gene,
    counts,
    check.names = FALSE
  ),
  file.path(result_dir, "counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  sample_data,
  file.path(result_dir, "sample_design.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

fit_start <- proc.time()
guide_result <- bb_screen(
  counts = counts,
  totals = colSums(counts),
  data = sample_data,
  formula = ~ marker_z + donor,
  term = "marker_z",
  guide = raw$sgRNA,
  gene = raw$Gene,
  min_total_count = 30,
  ncores = as.integer(Sys.getenv("BARCS_NCORES", "4"))
)
guide_seconds <- unname((proc.time() - fit_start)[["elapsed"]])
negative_control <- raw$Gene == "Non-Targeting Control"
calibrated <- bb_calibrate_controls(
  guide_result,
  negative_control,
  alpha = 0.05
)

aggregation_start <- proc.time()
original <- bb_gene_original(calibrated)
original_seconds <- unname((proc.time() - aggregation_start)[["elapsed"]])
aggregation_start <- proc.time()
normal <- bb_gene_normal(
  guide_result,
  min_guides = 3L,
  reference = "student_t"
)
normal_seconds <- unname((proc.time() - aggregation_start)[["elapsed"]])
aggregation_start <- proc.time()
partial <- bb_gene_partial_pool(
  guide_result,
  control = negative_control,
  min_guides = 2L,
  min_control_genes = 10L
)
partial_seconds <- unname((proc.time() - aggregation_start)[["elapsed"]])
aggregation_start <- proc.time()
eb <- bb_gene_eb_moderate(
  guide_result,
  control = negative_control,
  min_guides = 2L,
  prior_df = 4,
  min_control_genes = 10L
)
eb_seconds <- unname((proc.time() - aggregation_start)[["elapsed"]])

method_results <- list(
  `BARCS-original` = original,
  `BARCS-NORM` = normal,
  `BARCS-partial` = partial,
  `BARCS-EB` = eb
)
method_colours <- c(
  `BARCS-original` = "#0072B2",
  `BARCS-NORM` = "#7A3E9D",
  `BARCS-partial` = "#009E73",
  `BARCS-EB` = "#D55E00"
)
write.csv(
  guide_result,
  gzfile(file.path(result_dir, "guide_results.csv.gz")),
  row.names = FALSE
)
write.csv(
  calibrated,
  gzfile(file.path(result_dir, "calibrated_guide_results.csv.gz")),
  row.names = FALSE
)
for (method in names(method_results)) {
  write.csv(
    method_results[[method]],
    file.path(
      result_dir,
      paste0(gsub("-", "_", tolower(method)), "_gene_results.csv")
    ),
    row.names = FALSE
  )
}
write.csv(
  data.frame(
    analysis_protocol = analysis_protocol,
    method = names(method_results),
    shared_guide_fit_seconds = guide_seconds,
    gene_aggregation_seconds = c(
      original_seconds, normal_seconds, partial_seconds, eb_seconds
    ),
    guide_control_scale = attr(calibrated, "control_scale"),
    gene_null_scale = c(
      NA_real_, NA_real_,
      attr(partial, "null_scale"), attr(eb, "null_scale")
    ),
    prior_tau2 = c(
      NA_real_, NA_real_, NA_real_, attr(eb, "prior_tau2")
    ),
    prior_df = c(
      NA_real_, NA_real_, NA_real_, attr(eb, "prior_df")
    )
  ),
  file.path(result_dir, "runtime.csv"),
  row.names = FALSE
)

validation_raw <- read.csv(validation_source, check.names = FALSE)
validation_panel <- validation_raw[
  validation_raw$Gene_knocked_out != "Non-Targeting" &
    !is.na(validation_raw$Validation_status),
  ,
  drop = FALSE
]
validation_panel$gene <- validation_panel$Gene_knocked_out
validation_panel$truth <- grepl(
  "^Concordant", validation_panel$Validation_status
)
validation_panel$expected_sign <- ifelse(
  grepl("Increase", validation_panel$Validation_status),
  1,
  ifelse(grepl("Decrease", validation_panel$Validation_status), -1, NA)
)
if (sum(validation_panel$truth) != 26L ||
    sum(!validation_panel$truth) != 7L) {
  stop("Expected 26 validated and seven non-validating follow-up genes.")
}

safe_ratio <- function(numerator, denominator) {
  if (denominator > 0) numerator / denominator else NA_real_
}

validation_auc <- function(truth, score) {
  score_rank <- rank(score, ties.method = "average")
  n_positive <- sum(truth)
  n_negative <- sum(!truth)
  (
    sum(score_rank[truth]) -
      n_positive * (n_positive + 1) / 2
  ) / (n_positive * n_negative)
}

validation_average_precision <- function(truth, score) {
  ordered_truth <- truth[order(score, decreasing = TRUE)]
  precision_at_rank <- cumsum(ordered_truth) / seq_along(ordered_truth)
  mean(precision_at_rank[ordered_truth])
}

evaluate <- function(method, gene_result) {
  screen <- gene_result[
    gene_result$gene != "Non-Targeting Control" &
      is.finite(gene_result$p_value) &
      is.finite(gene_result$fdr),
    ,
    drop = FALSE
  ]
  assessed <- merge(
    validation_panel[
      , c("gene", "truth", "expected_sign", "Validation_status")
    ],
    screen[, c("gene", "estimate", "p_value", "fdr")],
    by = "gene",
    all.x = TRUE
  )
  assessed$called <- is.finite(assessed$fdr) & assessed$fdr < 0.10
  assessed$sign_match <- sign(assessed$estimate) == assessed$expected_sign
  assessed$directional_call <-
    assessed$called & assessed$sign_match & assessed$truth
  assessed$rank <- rank(
    screen$p_value, ties.method = "average"
  )[match(assessed$gene, screen$gene)]
  assessed$rank_percentile <- assessed$rank / nrow(screen)
  assessed$method <- method
  write.csv(
    assessed,
    file.path(
      result_dir,
      paste0(gsub("-", "_", tolower(method)), "_validation_panel.csv")
    ),
    row.names = FALSE
  )

  true_positive <- sum(assessed$called & assessed$truth)
  false_positive <- sum(assessed$called & !assessed$truth)
  true_negative <- sum(!assessed$called & !assessed$truth)
  false_negative <- sum(!assessed$called & assessed$truth)
  precision <- safe_ratio(
    true_positive, true_positive + false_positive
  )
  recall <- safe_ratio(
    true_positive, true_positive + false_negative
  )
  specificity <- safe_ratio(
    true_negative, true_negative + false_positive
  )
  f1 <- if (is.finite(precision) && is.finite(recall)) {
    safe_ratio(2 * precision * recall, precision + recall)
  } else {
    NA_real_
  }
  score <- -log10(pmax(assessed$p_value, .Machine$double.xmin))
  data.frame(
    analysis_protocol = analysis_protocol,
    method = method,
    screen_genes = nrow(screen),
    discoveries_at_fdr_0_10 = sum(screen$fdr < 0.10),
    directionally_validated_recovered =
      sum(assessed$directional_call, na.rm = TRUE),
    validated_total = sum(assessed$truth),
    validation_panel_true_positive = true_positive,
    validation_panel_false_positive = false_positive,
    validation_panel_true_negative = true_negative,
    validation_panel_false_negative = false_negative,
    precision = precision,
    recall = recall,
    specificity = specificity,
    balanced_accuracy = mean(c(recall, specificity), na.rm = TRUE),
    f1 = f1,
    auroc = validation_auc(assessed$truth, score),
    average_precision =
      validation_average_precision(assessed$truth, score),
    direction_concordance_validated =
      mean(assessed$sign_match[assessed$truth], na.rm = TRUE),
    median_validated_rank_percentile = median(
      assessed$rank_percentile[assessed$truth],
      na.rm = TRUE
    )
  )
}

metrics <- do.call(rbind, lapply(names(method_results), function(method) {
  evaluate(method, method_results[[method]])
}))
write.csv(
  metrics,
  file.path(result_dir, "benchmark_metrics.csv"),
  row.names = FALSE
)
write.csv(
  metrics,
  file.path(
    "data", "derived", "waterbear_facs_three_method_metrics.csv"
  ),
  row.names = FALSE
)

null_p <- guide_result$p_value[
  negative_control & is.finite(guide_result$p_value)
]
calibrated_null_p <- calibrated$p_value[
  negative_control & is.finite(calibrated$p_value)
]
null_metrics <- data.frame(
  analysis_protocol = analysis_protocol,
  guide_statistic = c("raw", "control_calibrated"),
  n_non_targeting_guides = c(length(null_p), length(calibrated_null_p)),
  fraction_p_below_0_05 = c(
    mean(null_p < 0.05),
    mean(calibrated_null_p < 0.05)
  ),
  fraction_p_below_0_10 = c(
    mean(null_p < 0.10),
    mean(calibrated_null_p < 0.10)
  ),
  median_p_value = c(median(null_p), median(calibrated_null_p))
)
write.csv(
  null_metrics,
  file.path(result_dir, "non_targeting_guide_calibration.csv"),
  row.names = FALSE
)
write.csv(
  null_metrics,
  file.path(
    "data", "derived",
    "waterbear_facs_three_method_null_calibration.csv"
  ),
  row.names = FALSE
)

pdf(figure_path, width = 10.5, height = 4.2, useDingbats = FALSE)
layout(matrix(1:3, nrow = 1))

par(mar = c(7, 4.2, 2.5, 0.8))
validation_positions <- barplot(
  metrics$directionally_validated_recovered,
  names.arg = metrics$method,
  col = unname(method_colours[metrics$method]),
  border = NA,
  las = 2,
  cex.names = 0.72,
  ylim = c(0, 26),
  ylab = "Genes recovered",
  main = "A  Directional validation"
)
abline(h = 26, lty = 3, col = "#555555")
text(
  validation_positions,
  metrics$directionally_validated_recovered,
  labels = metrics$directionally_validated_recovered,
  pos = 3,
  cex = 0.78
)

par(mar = c(7, 4.2, 2.5, 0.8))
panel_matrix <- rbind(
  F1 = metrics$f1,
  `Average precision` = metrics$average_precision,
  `Balanced accuracy` = metrics$balanced_accuracy
)
matplot(
  seq_len(nrow(panel_matrix)),
  panel_matrix,
  type = "b",
  pch = 16,
  lty = 1,
  lwd = 2,
  col = unname(method_colours[metrics$method]),
  xaxt = "n",
  xlab = "",
  las = 2,
  ylim = c(0, 1),
  ylab = "Selected-panel metric",
  main = "B  33-gene follow-up panel",
  bty = "l"
)
axis(
  1,
  at = seq_len(nrow(panel_matrix)),
  labels = rownames(panel_matrix),
  las = 2,
  cex.axis = 0.72
)
legend(
  "bottomleft",
  legend = metrics$method,
  col = unname(method_colours[metrics$method]),
  lty = 1,
  lwd = 2,
  pch = 16,
  bty = "n",
  cex = 0.68
)

par(mar = c(7, 4.2, 2.5, 0.8))
rank_upper <- max(
  0.03,
  1.2 * metrics$median_validated_rank_percentile
)
barplot(
  metrics$median_validated_rank_percentile,
  names.arg = metrics$method,
  col = unname(method_colours[metrics$method]),
  border = NA,
  las = 2,
  cex.names = 0.72,
  ylim = c(0, rank_upper),
  ylab = "Median rank / screen size",
  main = "C  Validated-gene rank\n(lower is better)"
)
dev.off()

print(metrics)
print(null_metrics)
