#!/usr/bin/env Rscript

# GSE242880 head-to-head comparison. The three current BARCS results and two
# independently rerun MAGeCK results are evaluated against the same selected
# 33-gene follow-up panel. Waterbear and MAUDE are retained as literature
# recovery references only because complete per-gene outputs for the panel
# are not available here.

options(stringsAsFactors = FALSE)
analysis_protocol <- "waterbear-external-headtohead-v1"
result_root <- file.path("results", "waterbear_facs")
barcs_root <- file.path(result_root, "three_methods")
validation_path <- file.path(
  "data", "raw", "waterbear", "cd25_validated_targets_revised.csv"
)

result_paths <- c(
  `BARCS-original` = file.path(
    barcs_root, "barcs_original_gene_results.csv"
  ),
  `BARCS-partial` = file.path(
    barcs_root, "barcs_partial_gene_results.csv"
  ),
  `BARCS-EB` = file.path(
    barcs_root, "barcs_eb_gene_results.csv"
  ),
  `MAGeCK-MLE` = file.path(
    result_root, "mageck_mle_all_bins_gene_results.csv"
  ),
  `MAGeCK-test` = file.path(
    result_root, "mageck_test_outer_bins_gene_results.csv"
  )
)
required <- c(validation_path, result_paths)
if (!all(file.exists(required))) {
  stop(
    "Missing GSE242880 validation or method result file(s): ",
    paste(required[!file.exists(required)], collapse = ", "),
    call. = FALSE
  )
}

validation_raw <- read.csv(validation_path, check.names = FALSE)
validation <- validation_raw[
  validation_raw$Gene_knocked_out != "Non-Targeting" &
    !is.na(validation_raw$Validation_status),
  ,
  drop = FALSE
]
validation$gene <- validation$Gene_knocked_out
validation$truth <- grepl("^Concordant", validation$Validation_status)
validation$expected_sign <- ifelse(
  grepl("Increase", validation$Validation_status),
  1,
  ifelse(grepl("Decrease", validation$Validation_status), -1, NA)
)
stopifnot(
  nrow(validation) == 33L,
  sum(validation$truth) == 26L,
  sum(!validation$truth) == 7L
)

safe_ratio <- function(numerator, denominator) {
  if (denominator > 0) numerator / denominator else NA_real_
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

matthews_correlation <- function(tp, fp, tn, fn) {
  denominator <- sqrt(
    (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
  )
  if (denominator > 0) {
    (tp * tn - fp * fn) / denominator
  } else {
    NA_real_
  }
}

evaluate <- function(method, gene_result, bins_used, design) {
  normalized_gene <- tolower(gsub("_", " ", gene_result$gene))
  screen <- gene_result[
    normalized_gene != "non-targeting control" &
      is.finite(gene_result$estimate) &
      is.finite(gene_result$p_value) &
      is.finite(gene_result$fdr),
    ,
    drop = FALSE
  ]
  assessed <- merge(
    validation[, c(
      "gene", "truth", "expected_sign", "Validation_status"
    )],
    screen[, c("gene", "estimate", "p_value", "fdr")],
    by = "gene",
    all.x = TRUE
  )
  if (any(!is.finite(assessed$p_value)) ||
      any(!is.finite(assessed$fdr)) ||
      any(!is.finite(assessed$estimate))) {
    stop(method, " is missing one or more validation-panel genes.")
  }
  assessed$called <- assessed$fdr < 0.10
  assessed$sign_match <-
    sign(assessed$estimate) == assessed$expected_sign
  assessed$directional_call <-
    assessed$called & assessed$sign_match & assessed$truth
  assessed$method <- method
  assessed$analysis_protocol <- analysis_protocol

  tp <- sum(assessed$called & assessed$truth)
  fp <- sum(assessed$called & !assessed$truth)
  tn <- sum(!assessed$called & !assessed$truth)
  fn <- sum(!assessed$called & assessed$truth)
  precision <- safe_ratio(tp, tp + fp)
  recall <- safe_ratio(tp, tp + fn)
  specificity <- safe_ratio(tn, tn + fp)
  f1 <- safe_ratio(2 * precision * recall, precision + recall)
  score <- -log10(pmax(assessed$p_value, .Machine$double.xmin))

  list(
    metric = data.frame(
      analysis_protocol = analysis_protocol,
      method = method,
      result_type = "independent rerun",
      bins_used = bins_used,
      design = design,
      screen_genes = nrow(screen),
      discoveries_at_fdr_0_10 = sum(screen$fdr < 0.10),
      directionally_validated_recovered =
        sum(assessed$directional_call),
      validated_total = sum(assessed$truth),
      validation_panel_true_positive = tp,
      validation_panel_false_positive = fp,
      validation_panel_true_negative = tn,
      validation_panel_false_negative = fn,
      precision = precision,
      recall = recall,
      specificity = specificity,
      balanced_accuracy = mean(c(recall, specificity)),
      matthews_correlation = matthews_correlation(tp, fp, tn, fn),
      f1 = f1,
      auroc = auroc(assessed$truth, score),
      average_precision =
        average_precision(assessed$truth, score),
      direction_concordance_validated =
        mean(assessed$sign_match[assessed$truth]),
      stringsAsFactors = FALSE
    ),
    panel = assessed
  )
}

method_design <- data.frame(
  method = names(result_paths),
  bins_used = c(4L, 4L, 4L, 4L, 2L),
  design = c(
    rep("ordered four-bin trend + donor", 4L),
    "outer-bin Q1 versus Q4"
  ),
  stringsAsFactors = FALSE
)

evaluated <- lapply(seq_len(nrow(method_design)), function(index) {
  method <- method_design$method[index]
  evaluate(
    method,
    read.csv(result_paths[[method]]),
    method_design$bins_used[index],
    method_design$design[index]
  )
})
rerun_metrics <- do.call(rbind, lapply(evaluated, `[[`, "metric"))
panel_results <- do.call(rbind, lapply(evaluated, `[[`, "panel"))

published_metrics <- data.frame(
  analysis_protocol = analysis_protocol,
  method = c("Waterbear", "MAUDE"),
  result_type = "published aggregate",
  bins_used = 4L,
  design = "joint four-bin FACS model",
  screen_genes = 1350L,
  discoveries_at_fdr_0_10 = c(79L, 406L),
  directionally_validated_recovered = c(24L, 25L),
  validated_total = 26L,
  validation_panel_true_positive = NA_integer_,
  validation_panel_false_positive = NA_integer_,
  validation_panel_true_negative = NA_integer_,
  validation_panel_false_negative = NA_integer_,
  precision = NA_real_,
  recall = NA_real_,
  specificity = NA_real_,
  balanced_accuracy = NA_real_,
  matthews_correlation = NA_real_,
  f1 = NA_real_,
  auroc = NA_real_,
  average_precision = NA_real_,
  direction_concordance_validated = NA_real_,
  stringsAsFactors = FALSE
)
metrics <- rbind(rerun_metrics, published_metrics)

write.csv(
  metrics,
  file.path(
    "data", "derived",
    "waterbear_facs_external_head_to_head_metrics.csv"
  ),
  row.names = FALSE
)
write.csv(
  panel_results,
  file.path(
    "data", "derived",
    "waterbear_facs_external_head_to_head_validation_panel.csv"
  ),
  row.names = FALSE
)

provenance <- data.frame(
  analysis_protocol = analysis_protocol,
  method = c(names(result_paths), "Waterbear", "MAUDE"),
  result_type = c(
    rep("independent rerun", length(result_paths)),
    "published aggregate", "published aggregate"
  ),
  source = c(
    unname(result_paths),
    "Pimentel et al. 2024, doi:10.1101/2024.06.17.599448",
    "Pimentel et al. 2024 comparison of MAUDE"
  ),
  md5 = c(
    unname(tools::md5sum(result_paths)),
    NA_character_, NA_character_
  ),
  per_gene_results_available = c(
    rep(TRUE, length(result_paths)), FALSE, FALSE
  ),
  stringsAsFactors = FALSE
)
write.csv(
  provenance,
  file.path(
    "data", "derived",
    "waterbear_facs_external_head_to_head_provenance.csv"
  ),
  row.names = FALSE
)

method_colours <- c(
  `BARCS-original` = "#0072B2",
  `BARCS-partial` = "#009E73",
  `BARCS-EB` = "#D55E00",
  `MAGeCK-MLE` = "#6A3D9A",
  `MAGeCK-test` = "#CC79A7",
  Waterbear = "#E69F00",
  MAUDE = "#999999"
)
pdf(
  file.path("figures", "waterbear_facs_external_head_to_head.pdf"),
  width = 11.5,
  height = 4.5,
  useDingbats = FALSE
)
layout(matrix(1:3, nrow = 1))

par(mar = c(7.6, 4.2, 2.7, 0.8))
positions <- barplot(
  metrics$directionally_validated_recovered,
  names.arg = metrics$method,
  col = unname(method_colours[metrics$method]),
  border = NA,
  las = 2,
  cex.names = 0.68,
  ylim = c(0, 28),
  ylab = "Genes recovered",
  main = "A  Directional recovery"
)
abline(h = 26, lty = 3, col = "#555555")
text(
  positions,
  metrics$directionally_validated_recovered,
  labels = paste0(
    metrics$directionally_validated_recovered, "/",
    metrics$validated_total
  ),
  pos = 3,
  cex = 0.72
)

par(mar = c(7.6, 4.2, 2.7, 0.8))
positions <- barplot(
  metrics$discoveries_at_fdr_0_10,
  names.arg = metrics$method,
  col = unname(method_colours[metrics$method]),
  border = NA,
  las = 2,
  cex.names = 0.68,
  ylim = c(0, 440),
  ylab = "Screen discoveries",
  main = "B  Calls at reported threshold"
)
text(
  positions,
  metrics$discoveries_at_fdr_0_10,
  labels = metrics$discoveries_at_fdr_0_10,
  pos = 3,
  cex = 0.72
)

par(mar = c(7.6, 4.2, 2.7, 0.8))
selected_metrics <- c("f1", "average_precision", "balanced_accuracy")
panel_matrix <- t(as.matrix(rerun_metrics[, selected_metrics]))
colnames(panel_matrix) <- rerun_metrics$method
barplot(
  panel_matrix,
  beside = TRUE,
  names.arg = rerun_metrics$method,
  col = c("#4C78A8", "#F2CF5B", "#59A14F"),
  border = NA,
  las = 2,
  cex.names = 0.68,
  ylim = c(0, 1.12),
  ylab = "Metric value",
  main = "C  Selected 33-gene panel"
)
legend(
  "top",
  legend = c("F1", "Average precision", "Balanced accuracy"),
  fill = c("#4C78A8", "#F2CF5B", "#59A14F"),
  border = NA,
  bty = "n",
  cex = 0.62,
  horiz = TRUE,
  inset = 0.01
)
dev.off()

print(metrics)
