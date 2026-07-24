#!/usr/bin/env Rscript

# Independent biological benchmark using the Sanson A375 screen bundled with
# CB2. Essential and nonessential reference genes define an observable target:
# a better ranking should separate the two classes more accurately.

source(file.path("R", "bbreg.R"))
source(file.path("R", "method_palette.R"))
load(file.path("CB2", "data", "Sanson_CRISPRn_A375.rda"))

benchmark_dir <- file.path("results", "sanson_benchmark")
dir.create(benchmark_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

counts <- as.matrix(Sanson_CRISPRn_A375$count)
storage.mode(counts) <- "double"
design <- Sanson_CRISPRn_A375$design
counts <- counts[, design$sample_name, drop = FALSE]
guide <- rownames(counts)
gene <- sub("_.*$", "", guide)
sample_data <- data.frame(after = as.integer(design$group == "trt"))

message("Fitting beta-binomial regression to ", nrow(counts), " guides ...")
bb_guide <- bb_screen(
  counts = counts,
  totals = colSums(counts),
  data = sample_data,
  formula = ~ after,
  term = "after",
  guide = guide,
  gene = gene,
  ncores = min(4L, parallel::detectCores(logical = FALSE))
)

gene_index <- split(seq_along(gene), gene)
bb_gene <- do.call(rbind, lapply(names(gene_index), function(gene_name) {
  index <- gene_index[[gene_name]]
  keep <- is.finite(bb_guide$p_value[index]) &
    is.finite(bb_guide$estimate[index])
  if (!any(keep)) {
    return(data.frame(
      gene = gene_name, estimate = NA_real_, p_value = NA_real_
    ))
  }
  index <- index[keep]
  signed_z <- sign(bb_guide$estimate[index]) *
    qnorm(
      pmax(bb_guide$p_value[index] / 2, .Machine$double.xmin),
      lower.tail = FALSE
    )
  combined_z <- sum(signed_z) / sqrt(length(signed_z))
  data.frame(
    gene = gene_name,
    estimate = median(bb_guide$estimate[index]),
    p_value = 2 * pnorm(-abs(combined_z))
  )
}))
bb_gene$fdr <- p.adjust(bb_gene$p_value, method = "BH")
bb_gene$score <- -sign(bb_gene$estimate) *
  -log10(pmax(bb_gene$p_value, .Machine$double.xmin))

cnv_path <- file.path("data", "derived", "A375_DepMap19Q3_CNV.tsv")
if (!file.exists(cnv_path)) {
  stop("A375 copy-number profile is missing at `", cnv_path, "`.", call. = FALSE)
}
cnv <- read.delim(cnv_path, stringsAsFactors = FALSE)
bb_effect_path <- file.path(benchmark_dir, "beta_binomial_gene_effects.tsv")
bb_cnv_effect_path <- file.path(
  benchmark_dir, "beta_binomial_gene_effects_cnv_corrected.tsv"
)
write.table(
  bb_gene[is.finite(bb_gene$estimate), c("gene", "estimate")],
  bb_effect_path,
  sep = "\t", quote = FALSE, row.names = FALSE
)

mageck_count_path <- file.path(benchmark_dir, "counts.tsv")
mageck_design_path <- file.path(benchmark_dir, "design.tsv")
mageck_prefix <- file.path(benchmark_dir, "mageck")
mageck_gene_path <- paste0(mageck_prefix, ".gene_summary.txt")
mageck_cnv_prefix <- file.path(benchmark_dir, "mageck_cnv")
mageck_cnv_gene_path <- paste0(mageck_cnv_prefix, ".gene_summary.txt")
write.table(
  data.frame(sgRNA = guide, Gene = gene, counts, check.names = FALSE),
  mageck_count_path,
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  data.frame(
    samples = design$sample_name,
    Initial_condition = 1,
    after = as.integer(design$group == "trt")
  ),
  mageck_design_path,
  sep = "\t", quote = FALSE, row.names = FALSE
)
mageck_executable <- file.path(".venv", "bin", "mageck")
mageck_python <- file.path(".venv", "bin", "python")
mageck_compat <- file.path("scripts", "mageck_compat.py")
mageck_cnv_correct <- file.path("scripts", "mageck_cnv_correct.py")
if (!all(file.exists(
  mageck_executable, mageck_python, mageck_compat, mageck_cnv_correct
))) {
  stop("Official MAGeCK is required at `.venv/bin/mageck`.", call. = FALSE)
}

bb_cnv_status <- system2(
  mageck_python,
  c(
    mageck_cnv_correct,
    "--effects", bb_effect_path,
    "--cnv", cnv_path,
    "--cell-line", "A375_SKIN",
    "--output", bb_cnv_effect_path
  )
)
if (bb_cnv_status != 0 || !file.exists(bb_cnv_effect_path)) {
  stop("MAGeCK piecewise CNV correction of beta-binomial effects failed.")
}
bb_cnv_effect <- read.delim(
  bb_cnv_effect_path, stringsAsFactors = FALSE
)
bb_gene$cnv_corrected_estimate <- bb_cnv_effect$cnv_corrected_estimate[
  match(bb_gene$gene, bb_cnv_effect$gene)
]
bb_gene$cnv_score <- -sign(bb_gene$cnv_corrected_estimate) *
  -log10(pmax(bb_gene$p_value, .Machine$double.xmin))

status <- system2(
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
if (status != 0 || !file.exists(mageck_gene_path)) {
  stop("Official MAGeCK test analysis failed.", call. = FALSE)
}
cnv_status <- system2(
  mageck_python,
  c(
    mageck_compat,
    "mle",
    "-k", mageck_count_path,
    "-d", mageck_design_path,
    "-n", mageck_cnv_prefix,
    "--norm-method", "median",
    "--permutation-round", "1",
    "--threads", "4",
    "--cnv-norm", cnv_path,
    "--cell-line", "A375_SKIN"
  )
)
if (cnv_status != 0 || !file.exists(mageck_cnv_gene_path)) {
  stop("Official MAGeCK CNV-corrected analysis failed.", call. = FALSE)
}
mageck_gene <- read.delim(
  mageck_gene_path, check.names = FALSE, stringsAsFactors = FALSE
)
mageck_cnv_gene <- read.delim(
  mageck_cnv_gene_path, check.names = FALSE, stringsAsFactors = FALSE
)

score <- merge(
  bb_gene[, c(
    "gene", "estimate", "cnv_corrected_estimate",
    "p_value", "fdr", "score", "cnv_score"
  )],
  data.frame(
    gene = mageck_gene$Gene,
    mageck_lfc = mageck_gene[["after|beta"]],
    mageck_p_value = mageck_gene[["after|wald-p-value"]],
    mageck_fdr = mageck_gene[["after|wald-fdr"]],
    mageck_score = -sign(mageck_gene[["after|beta"]]) *
      -log10(pmax(
        mageck_gene[["after|wald-p-value"]],
        .Machine$double.xmin
      ))
  ),
  by = "gene",
  all = FALSE
)
score <- merge(
  score,
  data.frame(
    gene = mageck_cnv_gene$Gene,
    mageck_cnv_lfc = mageck_cnv_gene[["after|beta"]],
    mageck_cnv_p_value = mageck_cnv_gene[["after|wald-p-value"]],
    mageck_cnv_fdr = mageck_cnv_gene[["after|wald-fdr"]],
    mageck_cnv_score = -sign(mageck_cnv_gene[["after|beta"]]) *
      -log10(pmax(
        mageck_cnv_gene[["after|wald-p-value"]],
        .Machine$double.xmin
      ))
  ),
  by = "gene",
  all = FALSE
)
score <- merge(
  score,
  cnv[, c("SYMBOL", "A375_SKIN")],
  by.x = "gene", by.y = "SYMBOL", all = FALSE
)
score$label <- ifelse(
  score$gene %in% Sanson_CRISPRn_A375$egenes, 1L,
  ifelse(score$gene %in% Sanson_CRISPRn_A375$ngenes, 0L, NA_integer_)
)
score <- score[!is.na(score$label), ]
nonessential <- score$label == 0L
cnv_bias <- data.frame(
  method = c(
    "beta-binomial regression",
    "beta-binomial + MAGeCK CNV correction",
    "official MAGeCK-MLE",
    "official MAGeCK-MLE + CNV"
  ),
  spearman_effect_vs_cnv_nonessential = c(
    cor(
      score$estimate[nonessential], score$A375_SKIN[nonessential],
      method = "spearman"
    ),
    cor(
      score$cnv_corrected_estimate[nonessential],
      score$A375_SKIN[nonessential],
      method = "spearman"
    ),
    cor(
      score$mageck_lfc[nonessential], score$A375_SKIN[nonessential],
      method = "spearman"
    ),
    cor(
      score$mageck_cnv_lfc[nonessential],
      score$A375_SKIN[nonessential],
      method = "spearman"
    )
  )
)

roc_points <- function(label, value) {
  order_index <- order(value, decreasing = TRUE, na.last = NA)
  y <- label[order_index]
  data.frame(
    fpr = c(0, cumsum(y == 0) / sum(y == 0)),
    tpr = c(0, cumsum(y == 1) / sum(y == 1))
  )
}
pr_points <- function(label, value) {
  order_index <- order(value, decreasing = TRUE, na.last = NA)
  y <- label[order_index]
  data.frame(
    recall = c(0, cumsum(y == 1) / sum(y == 1)),
    precision = c(1, cumsum(y == 1) / seq_along(y))
  )
}
auc_roc <- function(label, value) {
  keep <- is.finite(value)
  label <- label[keep]
  value <- value[keep]
  positive_rank_sum <- sum(rank(value, ties.method = "average")[label == 1])
  n_positive <- sum(label == 1)
  n_negative <- sum(label == 0)
  (positive_rank_sum - n_positive * (n_positive + 1) / 2) /
    (n_positive * n_negative)
}
average_precision <- function(label, value) {
  order_index <- order(value, decreasing = TRUE, na.last = NA)
  y <- label[order_index]
  precision <- cumsum(y == 1) / seq_along(y)
  mean(precision[y == 1])
}
recall_at_fpr <- function(label, value, target = 0.05) {
  curve <- roc_points(label, value)
  max(curve$tpr[curve$fpr <= target])
}
classification_metrics <- function(label, fdr, effect, threshold = 0.05) {
  predicted <- is.finite(fdr) & fdr <= threshold & effect < 0
  true_positive <- sum(predicted & label == 1)
  false_positive <- sum(predicted & label == 0)
  false_negative <- sum(!predicted & label == 1)
  precision <- if (true_positive + false_positive > 0) {
    true_positive / (true_positive + false_positive)
  } else {
    NA_real_
  }
  recall <- true_positive / (true_positive + false_negative)
  f1 <- if (is.finite(precision) && precision + recall > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    NA_real_
  }
  c(
    discoveries = true_positive + false_positive,
    precision = precision,
    recall = recall,
    F1 = f1
  )
}

method_score <- list(
  `beta-binomial regression` = score$score,
  `beta-binomial + MAGeCK CNV correction` = score$cnv_score,
  `official MAGeCK-MLE` = score$mageck_score,
  `official MAGeCK-MLE + CNV` = score$mageck_cnv_score
)
metrics <- do.call(rbind, lapply(names(method_score), function(method) {
  value <- method_score[[method]]
  data.frame(
    method = method,
    n_essential = sum(score$label == 1),
    n_nonessential = sum(score$label == 0),
    AUROC = auc_roc(score$label, value),
    average_precision = average_precision(score$label, value),
    recall_at_FPR_0.05 = recall_at_fpr(score$label, value)
  )
}))

fdr_threshold <- sort(unique(c(10^seq(-4, -1, length.out = 40L), 0.05)))
threshold_method <- list(
  `beta-binomial regression` = list(
    fdr = score$fdr, effect = score$estimate
  ),
  `beta-binomial + MAGeCK CNV correction` = list(
    fdr = score$fdr, effect = score$cnv_corrected_estimate
  ),
  `official MAGeCK-MLE` = list(
    fdr = score$mageck_fdr, effect = score$mageck_lfc
  ),
  `official MAGeCK-MLE + CNV` = list(
    fdr = score$mageck_cnv_fdr, effect = score$mageck_cnv_lfc
  )
)
threshold_metrics <- do.call(rbind, lapply(names(threshold_method), function(
  method
) {
  do.call(rbind, lapply(fdr_threshold, function(threshold) {
    value <- classification_metrics(
      score$label,
      threshold_method[[method]]$fdr,
      threshold_method[[method]]$effect,
      threshold
    )
    data.frame(
      method = method,
      fdr_threshold = threshold,
      discoveries = value[["discoveries"]],
      precision = value[["precision"]],
      recall = value[["recall"]],
      F1 = value[["F1"]]
    )
  }))
}))

set.seed(20260723)
n_bootstrap <- 2000L
positive <- which(score$label == 1)
negative <- which(score$label == 0)
bootstrap_difference <- replicate(n_bootstrap, {
  index <- c(
    sample(positive, length(positive), replace = TRUE),
    sample(negative, length(negative), replace = TRUE)
  )
  label <- score$label[index]
  bb_class <- classification_metrics(
    label, score$fdr[index], score$cnv_corrected_estimate[index], 0.05
  )
  mageck_class <- classification_metrics(
    label, score$mageck_cnv_fdr[index], score$mageck_cnv_lfc[index], 0.05
  )
  c(
    AUROC = auc_roc(label, score$cnv_score[index]) -
      auc_roc(label, score$mageck_cnv_score[index]),
    average_precision = average_precision(label, score$cnv_score[index]) -
      average_precision(label, score$mageck_cnv_score[index]),
    recall_at_FDR_0.05 = bb_class[["recall"]] -
      mageck_class[["recall"]],
    F1_at_FDR_0.05 = bb_class[["F1"]] - mageck_class[["F1"]]
  )
})
bb_class_005 <- classification_metrics(
  score$label, score$fdr, score$cnv_corrected_estimate, 0.05
)
mageck_class_005 <- classification_metrics(
  score$label, score$mageck_cnv_fdr, score$mageck_cnv_lfc, 0.05
)
bb_raw_class_005 <- classification_metrics(
  score$label, score$fdr, score$estimate, 0.05
)
mageck_raw_class_005 <- classification_metrics(
  score$label, score$mageck_fdr, score$mageck_lfc, 0.05
)
metric_value <- function(method, column) {
  metrics[metrics$method == method, column]
}
comparison <- data.frame(
  metric = c(
    "AUROC", "average_precision", "recall_at_FDR_0.05", "F1_at_FDR_0.05"
  ),
  beta_binomial_cnv_minus_mageck_cnv = c(
    metric_value(
      "beta-binomial + MAGeCK CNV correction", "AUROC"
    ) - metric_value("official MAGeCK-MLE + CNV", "AUROC"),
    metric_value(
      "beta-binomial + MAGeCK CNV correction", "average_precision"
    ) - metric_value("official MAGeCK-MLE + CNV", "average_precision"),
    bb_class_005[["recall"]] - mageck_class_005[["recall"]],
    bb_class_005[["F1"]] - mageck_class_005[["F1"]]
  ),
  bootstrap_lower_95 = apply(
    bootstrap_difference, 1, quantile, probs = 0.025, na.rm = TRUE
  ),
  bootstrap_upper_95 = apply(
    bootstrap_difference, 1, quantile, probs = 0.975, na.rm = TRUE
  )
)
cnv_change <- data.frame(
  method = c("beta-binomial regression", "official MAGeCK-MLE"),
  delta_AUROC = c(
    metric_value(
      "beta-binomial + MAGeCK CNV correction", "AUROC"
    ) - metric_value("beta-binomial regression", "AUROC"),
    metric_value(
      "official MAGeCK-MLE + CNV", "AUROC"
    ) - metric_value("official MAGeCK-MLE", "AUROC")
  ),
  delta_average_precision = c(
    metric_value(
      "beta-binomial + MAGeCK CNV correction", "average_precision"
    ) - metric_value("beta-binomial regression", "average_precision"),
    metric_value(
      "official MAGeCK-MLE + CNV", "average_precision"
    ) - metric_value("official MAGeCK-MLE", "average_precision")
  ),
  delta_recall_at_FDR_0.05 = c(
    bb_class_005[["recall"]] - bb_raw_class_005[["recall"]],
    mageck_class_005[["recall"]] - mageck_raw_class_005[["recall"]]
  ),
  delta_F1_at_FDR_0.05 = c(
    bb_class_005[["F1"]] - bb_raw_class_005[["F1"]],
    mageck_class_005[["F1"]] - mageck_raw_class_005[["F1"]]
  )
)

write.csv(
  score, file.path(benchmark_dir, "gold_standard_gene_scores.csv"),
  row.names = FALSE
)
write.csv(
  metrics, file.path(benchmark_dir, "benchmark_metrics.csv"),
  row.names = FALSE
)
write.csv(
  comparison, file.path(benchmark_dir, "paired_bootstrap_comparison.csv"),
  row.names = FALSE
)
write.csv(
  cnv_change, file.path(benchmark_dir, "cnv_adjustment_change.csv"),
  row.names = FALSE
)
write.csv(
  cnv_bias, file.path(benchmark_dir, "cnv_bias_diagnostic.csv"),
  row.names = FALSE
)
write.csv(
  threshold_metrics,
  file.path(benchmark_dir, "fdr_threshold_metrics.csv"),
  row.names = FALSE
)

roc_bb <- roc_points(score$label, score$score)
roc_bb_cnv <- roc_points(score$label, score$cnv_score)
roc_mageck <- roc_points(score$label, score$mageck_score)
roc_mageck_cnv <- roc_points(score$label, score$mageck_cnv_score)
pr_bb <- pr_points(score$label, score$score)
pr_bb_cnv <- pr_points(score$label, score$cnv_score)
pr_mageck <- pr_points(score$label, score$mageck_score)
pr_mageck_cnv <- pr_points(score$label, score$mageck_cnv_score)

pdf(file.path("figures", "sanson_gold_standard_benchmark.pdf"),
    width = 10, height = 8)
old_par <- par(mfrow = c(2, 2), mar = c(4.3, 4.4, 2.2, 1))
barcs_colour <- barcs_method_colours[["BARCS"]]
mageck_colour <- barcs_method_colours[["MAGeCK"]]
metric_labels <- c("BARCS", "BARCS + CNV", "MAGeCK-MLE", "MAGeCK-MLE + CNV")
plot(
  roc_bb$fpr, roc_bb$tpr, type = "l", lwd = 2, col = barcs_colour,
  xlab = "False-positive rate", ylab = "True-positive rate",
  main = "(A) Essential-gene ROC", xlim = c(0, 1), ylim = c(0, 1)
)
lines(roc_bb_cnv$fpr, roc_bb_cnv$tpr, lwd = 2, col = barcs_colour)
lines(roc_mageck$fpr, roc_mageck$tpr, lwd = 2, col = mageck_colour)
lines(
  roc_mageck_cnv$fpr, roc_mageck_cnv$tpr, lwd = 2, col = mageck_colour
)
lines(roc_bb$fpr, roc_bb$tpr, lwd = 1.5, lty = 2, col = barcs_colour)
lines(
  roc_mageck$fpr, roc_mageck$tpr, lwd = 1.5, lty = 2,
  col = mageck_colour
)
abline(0, 1, lty = 2, col = "grey50")
legend(
  "bottomright",
  legend = sprintf("%s (AUC %.3f)", metric_labels, metrics$AUROC),
  col = c(barcs_colour, barcs_colour, mageck_colour, mageck_colour),
  lty = c(2, 1, 2, 1), lwd = 2, bty = "n", cex = 0.55
)
plot(
  pr_bb$recall, pr_bb$precision, type = "l", lwd = 2,
  col = barcs_colour,
  xlab = "Recall", ylab = "Precision",
  main = "(B) Essential-gene precision-recall",
  xlim = c(0, 1), ylim = c(0, 1)
)
lines(pr_bb_cnv$recall, pr_bb_cnv$precision, lwd = 2, col = barcs_colour)
lines(pr_mageck$recall, pr_mageck$precision, lwd = 2, col = mageck_colour)
lines(
  pr_mageck_cnv$recall, pr_mageck_cnv$precision,
  lwd = 2, col = mageck_colour
)
lines(
  pr_bb$recall, pr_bb$precision, lwd = 1.5, lty = 2,
  col = barcs_colour
)
lines(
  pr_mageck$recall, pr_mageck$precision,
  lwd = 1.5, lty = 2, col = mageck_colour
)
abline(
  h = mean(score$label), lty = 2, col = "grey50"
)
legend(
  "bottomleft",
  legend = sprintf(
    "%s (AP %.3f)", metric_labels, metrics$average_precision
  ),
  col = c(barcs_colour, barcs_colour, mageck_colour, mageck_colour),
  lty = c(2, 1, 2, 1), lwd = 2, bty = "n", cex = 0.55
)
bb_threshold <- threshold_metrics[
  threshold_metrics$method == "beta-binomial regression", ]
bb_cnv_threshold <- threshold_metrics[
  threshold_metrics$method == "beta-binomial + MAGeCK CNV correction", ]
mageck_threshold <- threshold_metrics[
  threshold_metrics$method == "official MAGeCK-MLE", ]
mageck_cnv_threshold <- threshold_metrics[
  threshold_metrics$method == "official MAGeCK-MLE + CNV", ]
plot(
  bb_threshold$fdr_threshold, bb_threshold$F1,
  type = "l", log = "x", lwd = 2, col = barcs_colour,
  xlab = "Nominal FDR threshold", ylab = "Reference-gene F1",
  main = "(C) Thresholded discoveries",
  xlim = range(fdr_threshold), ylim = range(threshold_metrics$F1)
)
lines(
  bb_cnv_threshold$fdr_threshold, bb_cnv_threshold$F1,
  lwd = 2, col = barcs_colour
)
lines(
  mageck_threshold$fdr_threshold, mageck_threshold$F1,
  lwd = 1.5, lty = 2, col = mageck_colour
)
lines(
  mageck_cnv_threshold$fdr_threshold, mageck_cnv_threshold$F1,
  lwd = 2, col = mageck_colour
)
lines(
  bb_threshold$fdr_threshold, bb_threshold$F1,
  lwd = 1.5, lty = 2, col = barcs_colour
)
abline(v = 0.05, lty = 2, col = "grey50")
legend(
  "bottomright",
  legend = c(
    "BARCS", "BARCS + CNV",
    "official MAGeCK-MLE", "official MAGeCK-MLE + CNV"
  ),
  col = c(barcs_colour, barcs_colour, mageck_colour, mageck_colour),
  lty = c(2, 1, 2, 1), lwd = 2, bty = "n", cex = 0.58
)
cnv_values <- abs(cnv_bias$spearman_effect_vs_cnv_nonessential)
cnv_labels <- c("BARCS", "BARCS + CNV", "MAGeCK", "MAGeCK + CNV")
cnv_colours <- c(
  barcs_colour, barcs_colour, mageck_colour, mageck_colour
)
cnv_pch <- c(16, 1, 16, 1)
plot(
  0, 0, type = "n",
  xlim = c(0, max(cnv_values) * 1.08),
  ylim = c(0.5, length(cnv_values) + 0.5),
  xlab = "|Spearman(effect, copy number)| (closer to 0 is better)",
  ylab = "", yaxt = "n",
  main = "(D) Residual CNV association"
)
axis(
  2, at = seq_along(cnv_values), labels = rev(cnv_labels),
  las = 1, cex.axis = 0.78
)
abline(v = pretty(c(0, cnv_values)), col = "grey92", lwd = 0.8)
segments(
  x0 = 0, y0 = seq_along(cnv_values),
  x1 = rev(cnv_values), y1 = seq_along(cnv_values),
  col = rev(cnv_colours), lwd = 2
)
points(
  rev(cnv_values), seq_along(cnv_values),
  col = rev(cnv_colours), pch = rev(cnv_pch), cex = 1.15, lwd = 2
)
par(old_par)
dev.off()

cat("\nSanson gold-standard benchmark:\n")
print(metrics, row.names = FALSE, digits = 4)
cat("\nClassification at nominal FDR 0.05:\n")
print(
  rbind(
    `beta-binomial regression` = bb_raw_class_005,
    `beta-binomial + CNV` = bb_class_005,
    `official MAGeCK-MLE` = mageck_raw_class_005,
    `official MAGeCK-MLE + CNV` = mageck_class_005
  ),
  digits = 4
)
cat("\nPaired stratified-bootstrap difference:\n")
print(comparison, row.names = FALSE, digits = 4)
cat("\nChange after CNV correction:\n")
print(cnv_change, row.names = FALSE, digits = 4)
cat("\nResidual CNV association among reference nonessential genes:\n")
print(cnv_bias, row.names = FALSE, digits = 4)
