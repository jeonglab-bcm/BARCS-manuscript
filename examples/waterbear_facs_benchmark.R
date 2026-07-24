# Ordered-bin benchmark from the Waterbear IL2RA FACS screen (GSE242880).
#
# This is deliberately a scope-boundary benchmark.  The four FACS bins are
# ordered observations of a latent continuous IL2RA phenotype, but bins from
# one donor are a compositional partition rather than independent biological
# samples.  Waterbear models that joint structure; beta-binomial regression
# and MAGeCK-MLE do not.  We therefore report validated-hit recovery and
# calibration diagnostics without claiming that the regression is a
# replacement for a purpose-built FACS model.
#
# The low-coverage, high-MOI arm is the paper's most challenging real-data
# analysis.  It contains four bins for each of three donors and 26 genes whose
# direction was confirmed by individual knockout/flow-cytometry experiments.
# We compare:
#   1. beta-binomial regression using all four ordered bins;
#   2. official MAGeCK-MLE using the identical numeric-bin/donor design;
#   3. beta-binomial regression using only Q1 and Q4; and
#   4. official `mageck test` using Q1 versus Q4, as in the paper.
#
# Published Waterbear counts (24/26 validated genes; 79/1,350 discoveries) and
# MAUDE counts (25/26; 406/1,350) are included as labelled literature
# references.  Neither specialist method is rerun here.  Waterbear's deposited
# analysis uses four 10,000-iteration MCMC chains and the paper reports roughly
# 87 minutes for one fit.

options(stringsAsFactors = FALSE)

raw_dir <- file.path("data", "raw", "waterbear")
result_dir <- file.path("results", "waterbear_facs")
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
    "Download the GSE242880 processed count matrix and the Waterbear ",
    "validation table into `data/raw/waterbear` first.",
    call. = FALSE
  )
}

# Use the CB2 RcppArmadillo kernels when the local compiled library is present.
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

# The experiment targeted bin masses 0.2, 0.3, 0.3, and 0.2.  Conditional
# standard-normal means turn each interval into an interpretable latent-marker
# location instead of pretending that Q1, Q2, Q3, and Q4 are equally spaced.
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

combine_guides <- function(guide_result) {
  indices <- split(seq_len(nrow(guide_result)), guide_result$gene)
  result <- do.call(rbind, lapply(names(indices), function(gene_name) {
    index <- indices[[gene_name]]
    valid <- is.finite(guide_result$estimate[index]) &
      is.finite(guide_result$p_value[index])
    index <- index[valid]
    if (!length(index)) {
      return(NULL)
    }
    signed_z <- sign(guide_result$estimate[index]) *
      qnorm(
        pmax(guide_result$p_value[index] / 2, .Machine$double.xmin),
        lower.tail = FALSE
      )
    combined_z <- sum(signed_z) / sqrt(length(signed_z))
    data.frame(
      gene = gene_name,
      n_guides = length(index),
      estimate = median(guide_result$estimate[index]),
      p_value = 2 * pnorm(-abs(combined_z)),
      converged_fraction = mean(guide_result$converged[index])
    )
  }))
  result$fdr <- p.adjust(result$p_value, method = "BH")
  rownames(result) <- NULL
  result
}

run_bb <- function(count_matrix, design, cache_stem) {
  guide_path <- file.path(
    result_dir, paste0(cache_stem, "_guide_results.csv.gz")
  )
  if (!file.exists(guide_path) || identical(Sys.getenv("RERUN_BB"), "1")) {
    start <- proc.time()
    guide_result <- bb_screen(
      counts = count_matrix,
      totals = colSums(count_matrix),
      data = design,
      formula = ~ marker_z + donor,
      term = "marker_z",
      guide = raw$sgRNA,
      gene = raw$Gene,
      min_total_count = 30,
      ncores = 4
    )
    elapsed <- unname((proc.time() - start)[["elapsed"]])
    write.csv(guide_result, gzfile(guide_path), row.names = FALSE)
    write.csv(
      data.frame(
        method = cache_stem,
        n_guides = nrow(guide_result),
        elapsed_seconds = elapsed,
        guides_per_second = nrow(guide_result) / elapsed,
        converged = sum(guide_result$converged, na.rm = TRUE)
      ),
      file.path(result_dir, paste0(cache_stem, "_runtime.csv")),
      row.names = FALSE
    )
  } else {
    guide_result <- read.csv(gzfile(guide_path))
  }
  gene_result <- combine_guides(guide_result)
  write.csv(
    gene_result,
    file.path(result_dir, paste0(cache_stem, "_gene_results.csv")),
    row.names = FALSE
  )
  list(guide = guide_result, gene = gene_result)
}

bb_all <- run_bb(counts, sample_data, "beta_binomial_all_bins")
outer_index <- sample_data$bin %in% c(1L, 4L)
bb_outer <- run_bb(
  counts[, outer_index, drop = FALSE],
  droplevels(sample_data[outer_index, , drop = FALSE]),
  "beta_binomial_outer_bins"
)
negative_control <- raw$Gene == "Non-Targeting Control"
bb_all_calibrated_guide <- bb_calibrate_controls(
  bb_all$guide, negative_control, alpha = 0.05
)
bb_all_calibrated <- list(
  guide = bb_all_calibrated_guide,
  gene = combine_guides(bb_all_calibrated_guide)
)
write.csv(
  bb_all_calibrated$guide,
  gzfile(file.path(
    result_dir, "beta_binomial_all_bins_control_calibrated_guide_results.csv.gz"
  )),
  row.names = FALSE
)
write.csv(
  bb_all_calibrated$gene,
  file.path(
    result_dir, "beta_binomial_all_bins_control_calibrated_gene_results.csv"
  ),
  row.names = FALSE
)

mageck <- file.path(".venv", "bin", "mageck")
if (!file.exists(mageck)) {
  stop("Official MAGeCK 0.5.9.5 is required at `.venv/bin/mageck`.")
}
mageck_environment <- paste0(
  "PATH=", normalizePath(dirname(mageck)), ":", Sys.getenv("PATH")
)
# MAGeCK 0.5.9.5 tokenizes rows on all whitespace, so its input cannot retain
# the deposited "Non-Targeting Control" gene label verbatim.
mageck_count_path <- file.path(result_dir, "mageck_mle_counts.tsv")
mageck_keep <- rowSums(counts) > 0
write.table(
  data.frame(
    sgRNA = raw$sgRNA[mageck_keep],
    Gene = gsub("[[:space:]]+", "_", raw$Gene[mageck_keep]),
    counts[mageck_keep, , drop = FALSE],
    check.names = FALSE
  ),
  mageck_count_path,
  sep = "\t", quote = FALSE, row.names = FALSE
)
mageck_test_count_path <- file.path(result_dir, "mageck_test_counts.tsv")
write.table(
  data.frame(
    sgRNA = raw$sgRNA,
    Gene = gsub("[[:space:]]+", "_", raw$Gene),
    counts,
    check.names = FALSE
  ),
  mageck_test_count_path,
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Affine-equivalent all-bin continuous design used for beta-binomial
# regression.  MAGeCK-MLE's initializer requires at least one reference row
# whose non-intercept entries are all zero, so Q1 is mapped to 0 and Q4 to 1.
# An intercept makes this the same model space as the centered marker-z design.
mageck_design_path <- file.path(result_dir, "mageck_all_bin_design.tsv")
marker_range <- diff(range(sample_data$marker_z))
marker_unit <- (
  sample_data$marker_z - min(sample_data$marker_z)
) / marker_range
mageck_design <- data.frame(
  samples = sample_data$sample,
  baseline = 1,
  marker_unit = marker_unit,
  donor_D2 = as.integer(sample_data$donor == "D2"),
  donor_D3 = as.integer(sample_data$donor == "D3")
)
write.table(
  mageck_design,
  mageck_design_path,
  sep = "\t", quote = FALSE, row.names = FALSE
)
mageck_mle_prefix <- file.path(result_dir, "mageck_mle_all_bins")
mageck_mle_path <- paste0(mageck_mle_prefix, ".gene_summary.txt")
if (!file.exists(mageck_mle_path) ||
    identical(Sys.getenv("RERUN_MAGECK"), "1")) {
  status <- system2(
    mageck,
    c(
      "mle",
      "-k", mageck_count_path,
      "-d", mageck_design_path,
      "-n", mageck_mle_prefix,
      "--norm-method", "median",
      "--permutation-round", "1",
      "--no-permutation-by-group",
      "--threads", "4"
    ),
    env = mageck_environment
  )
  if (status != 0 || !file.exists(mageck_mle_path)) {
    stop("Official all-bin MAGeCK-MLE analysis failed.")
  }
}
mageck_mle_raw <- read.delim(mageck_mle_path, check.names = FALSE)
mageck_mle <- data.frame(
  gene = mageck_mle_raw$Gene,
  n_guides = mageck_mle_raw$sgRNA,
  estimate = mageck_mle_raw[["marker_unit|beta"]] / marker_range,
  p_value = mageck_mle_raw[["marker_unit|wald-p-value"]],
  fdr = mageck_mle_raw[["marker_unit|wald-fdr"]]
)
write.csv(
  mageck_mle,
  file.path(result_dir, "mageck_mle_all_bins_gene_results.csv"),
  row.names = FALSE
)

# Paper-matched MAGeCK analysis: low (Q1) treatment versus high (Q4) control.
q1 <- sample_data$sample[sample_data$bin == 1L]
q4 <- sample_data$sample[sample_data$bin == 4L]
mageck_test_prefix <- file.path(result_dir, "mageck_test_outer_bins")
mageck_test_path <- paste0(mageck_test_prefix, ".gene_summary.txt")
if (!file.exists(mageck_test_path) ||
    identical(Sys.getenv("RERUN_MAGECK"), "1")) {
  status <- system2(
    mageck,
    c(
      "test",
      "-k", mageck_test_count_path,
      "-t", paste(q1, collapse = ","),
      "-c", paste(q4, collapse = ","),
      "-n", mageck_test_prefix,
      "--norm-method", "median",
      "--sort-criteria", "pos"
    ),
    env = mageck_environment
  )
  if (status != 0 || !file.exists(mageck_test_path)) {
    stop("Official outer-bin MAGeCK test failed.")
  }
}
mageck_test_raw <- read.delim(mageck_test_path, check.names = FALSE)
mageck_test <- data.frame(
  gene = mageck_test_raw$id,
  n_guides = mageck_test_raw$num,
  # Treatment is Q1, so negate the reported low-versus-high fold change to
  # align positive effects with increasing IL2RA marker level.
  estimate = -mageck_test_raw[["pos|lfc"]],
  p_value = pmin(
    mageck_test_raw[["neg|p-value"]],
    mageck_test_raw[["pos|p-value"]]
  ),
  fdr = pmin(
    mageck_test_raw[["neg|fdr"]],
    mageck_test_raw[["pos|fdr"]]
  )
)
write.csv(
  mageck_test,
  file.path(result_dir, "mageck_test_outer_bins_gene_results.csv"),
  row.names = FALSE
)

validation_source_table <- read.csv(validation_source)
validation_panel <- validation_source_table[
  validation_source_table$Gene_knocked_out != "Non-Targeting" &
    !is.na(validation_source_table$Validation_status),
  ,
  drop = FALSE
]
validation_panel$gene <- validation_panel$Gene_knocked_out
validation_panel$truth <- grepl(
  "^Concordant", validation_panel$Validation_status
)
if (nrow(validation_panel) != 33L ||
    sum(validation_panel$truth) != 26L ||
    sum(!validation_panel$truth) != 7L) {
  stop("Expected a validation panel with 26 positives and 7 non-validating genes.")
}

validation <- validation_panel[validation_panel$truth, , drop = FALSE]
validation$gene <- validation$Gene_knocked_out
validation$expected_sign <- ifelse(
  grepl("Increase", validation$Validation_status), 1, -1
)
if (nrow(validation) != 26L) {
  stop("Expected the 26 directionally validated IL2RA regulators.")
}

evaluate_method <- function(method, gene_result) {
  gene_result <- gene_result[
    gene_result$gene != "Non-Targeting Control" &
      is.finite(gene_result$p_value) &
      is.finite(gene_result$fdr),
    ,
    drop = FALSE
  ]
  assessed <- merge(
    validation[, c("gene", "expected_sign")],
    gene_result[, c("gene", "estimate", "p_value", "fdr")],
    by = "gene",
    all.x = TRUE
  )
  assessed$sign_match <- sign(assessed$estimate) == assessed$expected_sign
  assessed$detected <- assessed$fdr < 0.10 & assessed$sign_match
  assessed$method <- method
  assessed$absolute_rank <- rank(
    gene_result$p_value, ties.method = "average"
  )[match(assessed$gene, gene_result$gene)]
  assessed$rank_percentile <- assessed$absolute_rank / nrow(gene_result)
  write.csv(
    assessed,
    file.path(
      result_dir,
      paste0(gsub("[^A-Za-z0-9]+", "_", tolower(method)),
             "_validated_genes.csv")
    ),
    row.names = FALSE
  )
  data.frame(
    method = method,
    bins_used = if (grepl("all", method, ignore.case = TRUE)) 4 else 2,
    independently_rerun = TRUE,
    discoveries_at_fdr_0_10 = sum(gene_result$fdr < 0.10),
    validated_recovered = sum(assessed$detected, na.rm = TRUE),
    validated_total = nrow(validation),
    validated_recall = mean(assessed$detected, na.rm = TRUE),
    direction_concordance = mean(assessed$sign_match, na.rm = TRUE),
    median_validated_rank_percentile =
      median(assessed$rank_percentile, na.rm = TRUE)
  )
}

benchmark <- rbind(
  evaluate_method(
    "Beta-binomial all bins + control calibration",
    bb_all_calibrated$gene
  ),
  evaluate_method("Beta-binomial all bins", bb_all$gene),
  evaluate_method("MAGeCK-MLE all bins", mageck_mle),
  evaluate_method("Beta-binomial outer bins", bb_outer$gene),
  evaluate_method("MAGeCK test outer bins", mageck_test),
  data.frame(
    method = "Waterbear published",
    bins_used = 4,
    independently_rerun = FALSE,
    discoveries_at_fdr_0_10 = 79,
    validated_recovered = 24,
    validated_total = 26,
    validated_recall = 24 / 26,
    direction_concordance = NA_real_,
    median_validated_rank_percentile = NA_real_
  ),
  data.frame(
    method = "MAUDE published",
    bins_used = 4,
    independently_rerun = FALSE,
    discoveries_at_fdr_0_10 = 406,
    validated_recovered = 25,
    validated_total = 26,
    validated_recall = 25 / 26,
    direction_concordance = NA_real_,
    median_validated_rank_percentile = NA_real_
  )
)
write.csv(
  benchmark,
  file.path(result_dir, "benchmark_metrics.csv"),
  row.names = FALSE
)

# The paper followed up 33 candidate genes individually: 26 had a concordant
# effect and seven did not validate.  This restricted panel provides actual
# positive and negative labels, unlike the rest of the screen, whose genes are
# unlabelled rather than false.  Metrics from this selected, small panel are
# useful supporting evidence but are not an unbiased genome-wide estimate.
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

evaluate_validation_panel <- function(method, gene_result) {
  assessed <- merge(
    validation_panel[, c("gene", "truth")],
    gene_result[, c("gene", "estimate", "p_value", "fdr")],
    by = "gene",
    all.x = TRUE
  )
  assessed$called <- is.finite(assessed$fdr) & assessed$fdr < 0.10
  tp <- sum(assessed$called & assessed$truth)
  fp <- sum(assessed$called & !assessed$truth)
  tn <- sum(!assessed$called & !assessed$truth)
  fn <- sum(!assessed$called & assessed$truth)
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  f1 <- 2 * precision * recall / (precision + recall)
  mcc <- (tp * tn - fp * fn) /
    sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  score <- -log10(pmax(assessed$p_value, .Machine$double.xmin))

  assessed$method <- method
  write.csv(
    assessed,
    file.path(
      result_dir,
      paste0(
        gsub("[^A-Za-z0-9]+", "_", tolower(method)),
        "_validation_panel.csv"
      )
    ),
    row.names = FALSE
  )

  data.frame(
    method = method,
    panel_positive = sum(assessed$truth),
    panel_negative = sum(!assessed$truth),
    true_positive = tp,
    false_positive = fp,
    true_negative = tn,
    false_negative = fn,
    precision = precision,
    recall = recall,
    specificity = specificity,
    balanced_accuracy = (recall + specificity) / 2,
    f1 = f1,
    mcc = mcc,
    auroc = validation_auc(assessed$truth, score),
    average_precision =
      validation_average_precision(assessed$truth, score)
  )
}

validation_metrics <- rbind(
  evaluate_validation_panel(
    "Beta-binomial all bins + control calibration",
    bb_all_calibrated$gene
  ),
  evaluate_validation_panel("Beta-binomial all bins", bb_all$gene),
  evaluate_validation_panel("MAGeCK-MLE all bins", mageck_mle),
  evaluate_validation_panel("Beta-binomial outer bins", bb_outer$gene),
  evaluate_validation_panel("MAGeCK test outer bins", mageck_test)
)
write.csv(
  validation_metrics,
  file.path(result_dir, "validation_panel_metrics.csv"),
  row.names = FALSE
)

# Non-targeting guides provide a direct null diagnostic for the two regression
# fits.  This is not a gene-level comparison with MAGeCK, whose RRA output does
# not report an equivalent guide-level p-value.
null_diagnostic <- function(method, guide_result) {
  p <- guide_result$p_value[
    guide_result$gene == "Non-Targeting Control" &
      is.finite(guide_result$p_value)
  ]
  data.frame(
    method = method,
    n_non_targeting_guides = length(p),
    fraction_p_below_0_05 = mean(p < 0.05),
    fraction_p_below_0_10 = mean(p < 0.10),
    median_p_value = median(p),
    genomic_inflation_lambda =
      median(qchisq(pmax(1 - p, 0), df = 1)) / qchisq(0.5, df = 1)
  )
}
null_metrics <- rbind(
  null_diagnostic(
    "Beta-binomial all bins + control calibration",
    bb_all_calibrated$guide
  ),
  null_diagnostic("Beta-binomial all bins", bb_all$guide),
  null_diagnostic("Beta-binomial outer bins", bb_outer$guide)
)
write.csv(
  null_metrics,
  file.path(result_dir, "non_targeting_calibration.csv"),
  row.names = FALSE
)

shared <- merge(
  bb_all$gene[, c("gene", "estimate", "p_value")],
  mageck_mle[, c("gene", "estimate", "p_value")],
  by = "gene",
  suffixes = c("_bb", "_mageck")
)
effect_correlation <- data.frame(
  comparison = "Beta-binomial all bins vs MAGeCK-MLE all bins",
  n_genes = nrow(shared),
  spearman_effect = cor(
    shared$estimate_bb, shared$estimate_mageck,
    method = "spearman", use = "complete.obs"
  ),
  spearman_significance = cor(
    -log10(shared$p_value_bb),
    -log10(shared$p_value_mageck),
    method = "spearman", use = "complete.obs"
  )
)
write.csv(
  effect_correlation,
  file.path(result_dir, "effect_rank_correlations.csv"),
  row.names = FALSE
)

method_colours <- c(
  "#006D5B", "#1B9E77", "#7570B3", "#66C2A5", "#8DA0CB", "#D95F02",
  "#E6AB02"
)
pdf(figure_path, width = 11, height = 8.5, useDingbats = FALSE)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
short_method <- c(
  "CB2-Reg 4-bin + ctrl",
  "CB2-Reg 4-bin raw",
  "MAGeCK-MLE 4-bin",
  "CB2-Reg outer-bin",
  "MAGeCK test outer",
  "Waterbear (reported)",
  "MAUDE (reported)"
)
par(mar = c(4, 10, 3, 1))
barplot(
  benchmark$validated_recovered,
  names.arg = short_method,
  horiz = TRUE,
  las = 1,
  col = method_colours,
  border = NA,
  xlim = c(0, 26),
  xlab = "Directionally validated genes recovered",
  main = "A. Positive-control recovery at 0.10 threshold"
)
abline(v = 26, lty = 3, col = "grey40")

panel_colours <- method_colours[seq_len(nrow(validation_metrics))]
barplot(
  validation_metrics$f1,
  names.arg = short_method[seq_len(nrow(validation_metrics))],
  horiz = TRUE,
  las = 1,
  col = panel_colours,
  border = NA,
  xlim = c(0, 1),
  xlab = "F1 in selected 33-gene validation panel",
  main = "B. Restricted validation-panel classification"
)

par(mar = c(5, 5, 3, 1))
plot(
  shared$estimate_mageck,
  shared$estimate_bb,
  pch = 16,
  cex = 0.45,
  col = adjustcolor("grey25", alpha.f = 0.25),
  xlab = "MAGeCK-MLE marker-z effect",
  ylab = "CB2-Reg marker-z effect",
  main = sprintf(
    "C. All-bin effects (Spearman rho = %.3f)",
    effect_correlation$spearman_effect
  )
)
abline(h = 0, v = 0, col = "grey75", lty = 3)
validated_index <- match(validation$gene, shared$gene)
validated_index <- validated_index[!is.na(validated_index)]
points(
  shared$estimate_mageck[validated_index],
  shared$estimate_bb[validated_index],
  pch = 21, bg = "#E7298A", col = "white", cex = 0.9
)

null_all <- sort(
  bb_all$guide$p_value[
    bb_all$guide$gene == "Non-Targeting Control" &
      is.finite(bb_all$guide$p_value)
  ]
)
null_all_calibrated <- sort(
  bb_all_calibrated$guide$p_value[
    bb_all_calibrated$guide$gene == "Non-Targeting Control" &
      is.finite(bb_all_calibrated$guide$p_value)
  ]
)
null_outer <- sort(
  bb_outer$guide$p_value[
    bb_outer$guide$gene == "Non-Targeting Control" &
      is.finite(bb_outer$guide$p_value)
  ]
)
plot(
  c(0, 1), c(0, 1),
  type = "n",
  xlab = "Nominal p-value",
  ylab = "Empirical cumulative fraction",
  main = "D. Non-targeting-guide null calibration"
)
abline(0, 1, lty = 3, col = "grey50")
lines(null_all, seq_along(null_all) / length(null_all),
      col = method_colours[2], lwd = 2)
lines(
  null_all_calibrated,
  seq_along(null_all_calibrated) / length(null_all_calibrated),
  col = method_colours[1], lwd = 2
)
lines(null_outer, seq_along(null_outer) / length(null_outer),
      col = method_colours[4], lwd = 2)
legend(
  "bottomright",
  legend = c(
    "CB2-Reg all bins + controls",
    "CB2-Reg all bins raw",
    "CB2-Reg outer bins raw"
  ),
  col = method_colours[c(1, 2, 4)],
  lwd = 2,
  bty = "n"
)
dev.off()

print(benchmark)
print(validation_metrics)
print(null_metrics)
print(effect_correlation)
