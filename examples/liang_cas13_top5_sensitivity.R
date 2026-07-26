#!/usr/bin/env Rscript

# Outcome-independent top-five-guide sensitivity analysis for Liang et al.
#
# The five guides are selected within each gene by their mean Day-0 abundance,
# before Day-14 depletion or any guide p-value is examined. Selecting guides by
# the smallest observed p-values would be circular and is intentionally not
# implemented. The primary Liang benchmark continues to use every valid guide.

options(stringsAsFactors = FALSE)
source(file.path("R", "bbreg.R"))

raw_dir <- file.path("data", "raw", "liang_cas13")
result_dir <- file.path("results", "liang_cas13")
derived_dir <- file.path("data", "derived")
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)

guide <- read.delim(
  file.path(raw_dir, "guide_library.tsv"),
  check.names = FALSE
)
expression <- read.delim(
  file.path(raw_dir, "lncrna_expression.tsv"),
  check.names = FALSE
)
cell_lines <- c("HAP1", "HEK293FT", "K562", "MDA-MB-231", "THP1")

roc_auc <- function(label, score) {
  rank_value <- rank(score, ties.method = "average")
  positives <- sum(label == 1L)
  negatives <- sum(label == 0L)
  (sum(rank_value[label == 1L]) -
     positives * (positives + 1) / 2) / (positives * negatives)
}

average_precision <- function(label, score) {
  label <- label[order(score, decreasing = TRUE)]
  precision <- cumsum(label == 1L) / seq_along(label)
  mean(precision[label == 1L])
}

recall_at_null_fpr <- function(label, score, target_fpr = 0.05) {
  cutoff <- unname(quantile(
    score[label == 0L],
    probs = 1 - target_fpr,
    type = 1
  ))
  mean(score[label == 1L] > cutoff)
}

gene_group <- unique(guide[, c("gene", "target_group")])
group_priority <- match(
  gene_group$target_group,
  c(
    "essential protein-coding gene", "long non-coding RNA",
    "protein-coding gene", "non-targeting"
  )
)
gene_group <- gene_group[order(group_priority), ]
gene_group <- gene_group[!duplicated(gene_group$gene), ]

evaluate <- function(gene_result, cell_line, rule) {
  gene_result <- merge(
    gene_result, gene_group,
    by = "gene", all.x = TRUE
  )
  gene_result <- merge(
    gene_result,
    expression[
      expression$cell_line == cell_line,
      c("gene", "tpm")
    ],
    by = "gene", all.x = TRUE
  )
  gene_result$truth <- ifelse(
    gene_result$target_group == "essential protein-coding gene", 1L,
    ifelse(
      gene_result$target_group == "long non-coding RNA" &
        is.finite(gene_result$tpm) & gene_result$tpm == 0,
      0L,
      NA_integer_
    )
  )
  keep <- !is.na(gene_result$truth) &
    is.finite(gene_result$estimate) &
    is.finite(gene_result$p_value)
  x <- gene_result[keep, , drop = FALSE]
  depletion_score <- -sign(x$estimate) *
    -log10(pmax(x$p_value, .Machine$double.xmin))
  data.frame(
    cell_line = cell_line,
    guide_rule = rule,
    positives = sum(x$truth == 1L),
    nulls = sum(x$truth == 0L),
    mean_guides_per_gene = mean(x$n_guides),
    auroc = roc_auc(x$truth, depletion_score),
    average_precision = average_precision(x$truth, depletion_score),
    recall_at_5pct_null_fpr = recall_at_null_fpr(
      x$truth, depletion_score
    ),
    null_p_lt_0_05 = mean(x$p_value[x$truth == 0L] < 0.05),
    essential_fdr_0_10_recall = mean(
      x$fdr[x$truth == 1L] < 0.10 &
        x$estimate[x$truth == 1L] < 0
    )
  )
}

metric_rows <- list()
row_index <- 0L
for (cell_line in cell_lines) {
  stem <- gsub("-", "_", cell_line)
  guide_path <- file.path(
    result_dir, paste0(stem, "_barcs_guide.csv.gz")
  )
  if (!file.exists(guide_path)) {
    stop(
      "Run `examples/liang_cas13_benchmark.R` before this sensitivity."
    )
  }
  guide_result <- read.csv(gzfile(guide_path))
  guide_index <- match(guide_result$guide, guide$sgrna)
  cell_guide <- guide[guide_index, , drop = FALSE]
  guide_result <- bb_calibrate_controls(
    guide_result,
    cell_guide$target_group == "non-targeting"
  )

  processed <- read.delim(
    file.path(
      raw_dir,
      paste0("published_processed_counts_", stem, ".tsv")
    ),
    check.names = FALSE
  )
  day0_columns <- grep(
    "^Day0 Replicate [12] +\\(Count\\)$",
    names(processed),
    value = TRUE
  )
  baseline <- rowMeans(processed[, day0_columns, drop = FALSE])
  names(baseline) <- processed$sgrna
  guide_result$baseline_abundance <- baseline[guide_result$guide]

  testable <- guide_result$gene != "non-targeting"
  all_gene <- bb_gene_original(
    guide_result[testable, , drop = FALSE],
    min_guides = 1L
  )

  groups <- split(
    which(testable),
    guide_result$gene[testable]
  )
  top5_index <- unlist(lapply(groups, function(index) {
    ordered <- index[order(
      guide_result$baseline_abundance[index],
      decreasing = TRUE,
      na.last = TRUE
    )]
    ordered[seq_len(min(5L, length(ordered)))]
  }), use.names = FALSE)
  top5_gene <- bb_gene_original(
    guide_result[top5_index, , drop = FALSE],
    min_guides = 1L
  )

  row_index <- row_index + 1L
  metric_rows[[row_index]] <- evaluate(
    all_gene, cell_line, "all_valid_guides"
  )
  row_index <- row_index + 1L
  metric_rows[[row_index]] <- evaluate(
    top5_gene, cell_line, "top5_by_mean_day0_abundance"
  )
}

metrics <- do.call(rbind, metric_rows)
macro <- aggregate(
  metrics[, c(
    "mean_guides_per_gene", "auroc", "average_precision",
    "recall_at_5pct_null_fpr", "null_p_lt_0_05",
    "essential_fdr_0_10_recall"
  )],
  by = list(guide_rule = metrics$guide_rule),
  FUN = mean
)

write.csv(
  metrics,
  file.path(derived_dir, "liang_cas13_top5_sensitivity_by_cell_line.csv"),
  row.names = FALSE
)
write.csv(
  macro,
  file.path(derived_dir, "liang_cas13_top5_sensitivity_macro_average.csv"),
  row.names = FALSE
)

print(metrics, row.names = FALSE)
cat("\nMacro-average:\n")
print(macro, row.names = FALSE)
