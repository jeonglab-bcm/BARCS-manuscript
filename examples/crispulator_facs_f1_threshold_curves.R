#!/usr/bin/env Rscript

# F1 and calibration curves across nominal gene-FDR thresholds for the eight
# methods in the CRISPulator low + bulk + high comparison.

options(stringsAsFactors = FALSE)
analysis_protocol <- "crispulator-f1-threshold-curves-v1"
thresholds <- c(0.10, 0.05, 0.01, 0.005, 0.001)
seeds <- 20250724L:20250728L
scenario_table <- read.csv(
  file.path("data", "derived", "crispulator_facs_parameter_grid.csv")
)
methods <- c(
  "BARCS-original", "BARCS-NORM", "BARCS-partial", "BARCS-EB",
  "MAGeCK-MLE", "edgeR-QL", "DESeq2", "limma-voom"
)
method_file <- c(
  `BARCS-original` = "low_bulk_high_original_gene_results.csv",
  `BARCS-NORM` = "low_bulk_high_normal_gene_results.csv",
  `BARCS-partial` = "low_bulk_high_partial_gene_results.csv",
  `BARCS-EB` = "low_bulk_high_eb_gene_results.csv",
  `MAGeCK-MLE` = "mageck_mle_low_bulk_high_gene_results.csv",
  `edgeR-QL` = "edger_ql_low_bulk_high_gene_results.csv",
  DESeq2 = "deseq2_low_bulk_high_gene_results.csv",
  `limma-voom` = "limma_voom_low_bulk_high_gene_results.csv"
)

format_parameter <- function(value) {
  format(value, trim = TRUE, scientific = FALSE)
}

cache_scenario_name <- function(moi, quality, genes, replicates) {
  name <- sprintf(
    "moi_%s_quality_%s",
    format_parameter(moi),
    format_parameter(quality)
  )
  if (genes != 400L) {
    name <- paste0(name, "_genes_", genes)
  }
  if (replicates != 4L) {
    name <- paste0(name, "_replicates_", replicates)
  }
  name
}

safe_ratio <- function(numerator, denominator, empty = NA_real_) {
  if (denominator > 0) numerator / denominator else empty
}

evaluate_threshold <- function(
  method,
  gene_result,
  gene_truth,
  threshold
) {
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
  active <- assessed$active
  called <- assessed$fdr < threshold
  true_positive <- sum(called & active)
  false_positive <- sum(called & !active)
  true_negative <- sum(!called & !active)
  false_negative <- sum(!called & active)
  precision <- safe_ratio(
    true_positive,
    true_positive + false_positive,
    empty = 0
  )
  recall <- safe_ratio(
    true_positive,
    true_positive + false_negative,
    empty = 0
  )
  f1 <- safe_ratio(
    2 * precision * recall,
    precision + recall,
    empty = 0
  )
  sign_match <- sign(assessed$estimate) == assessed$expected_sign
  data.frame(
    analysis_protocol = analysis_protocol,
    method = method,
    nominal_fdr = threshold,
    genes = nrow(assessed),
    active_genes = sum(active),
    calls = sum(called),
    true_positive = true_positive,
    false_positive = false_positive,
    true_negative = true_negative,
    false_negative = false_negative,
    precision = precision,
    recall = recall,
    directional_recall = mean(called[active] & sign_match[active]),
    f1 = f1,
    realized_fdp = safe_ratio(
      false_positive,
      true_positive + false_positive,
      empty = 0
    ),
    specificity = safe_ratio(
      true_negative,
      true_negative + false_positive
    ),
    balanced_accuracy = mean(c(
      recall,
      safe_ratio(true_negative, true_negative + false_positive)
    ))
  )
}

rows <- list()
row_index <- 0L
for (scenario_index in seq_len(nrow(scenario_table))) {
  scenario <- scenario_table[scenario_index, ]
  cache_scenario <- cache_scenario_name(
    scenario$moi,
    scenario$high_quality_guide_fraction,
    scenario$genes,
    scenario$replicates
  )
  for (seed in seeds) {
    run_root <- file.path(
      "results", "crispulator_facs", "repeated",
      cache_scenario, paste0("seed_", seed)
    )
    truth_path <- file.path(run_root, "data", "gene_truth.tsv")
    current_dir <- file.path(run_root, "analysis_three_methods")
    external_dir <- file.path(run_root, "analysis")
    result_paths <- c(
      file.path(
        current_dir,
        method_file[c(
          "BARCS-original", "BARCS-NORM", "BARCS-partial", "BARCS-EB"
        )]
      ),
      file.path(
        external_dir,
        method_file[c(
          "MAGeCK-MLE", "edgeR-QL", "DESeq2", "limma-voom"
        )]
      )
    )
    names(result_paths) <- methods
    required <- c(truth_path, result_paths)
    if (!all(file.exists(required))) {
      stop(
        "Missing truth or gene-result file(s) for ",
        scenario$scenario_id, ", seed ", seed, "."
      )
    }

    gene_truth <- read.delim(truth_path)
    gene_truth$active <-
      tolower(as.character(gene_truth$active)) == "true"
    method_results <- lapply(result_paths, read.csv)
    for (method in methods) {
      for (threshold in thresholds) {
        row_index <- row_index + 1L
        metric <- evaluate_threshold(
          method,
          method_results[[method]],
          gene_truth,
          threshold
        )
        metric$scenario_id <- scenario$scenario_id
        metric$seed <- seed
        metric$moi <- scenario$moi
        metric$high_quality_guide_fraction <-
          scenario$high_quality_guide_fraction
        metric$genes_simulated <- scenario$genes
        metric$replicates <- scenario$replicates
        metric$is_baseline <- scenario$is_baseline
        rows[[row_index]] <- metric
      }
    }
  }
}
metrics <- do.call(rbind, rows)
stopifnot(
  nrow(metrics) ==
    nrow(scenario_table) * length(seeds) * length(methods) *
      length(thresholds),
  setequal(unique(metrics$method), methods),
  setequal(unique(metrics$nominal_fdr), thresholds),
  all(metrics$f1 >= 0 & metrics$f1 <= 1),
  all(metrics$realized_fdp >= 0 & metrics$realized_fdp <= 1)
)
write.csv(
  metrics,
  file.path(
    "data", "derived", "crispulator_facs_f1_by_fdr.csv"
  ),
  row.names = FALSE
)

summarize_metrics <- function(data, grouping) {
  metric_names <- c(
    "precision", "recall", "directional_recall", "f1",
    "realized_fdp", "specificity", "balanced_accuracy", "calls"
  )
  group_id <- interaction(
    data[, grouping, drop = FALSE],
    drop = TRUE,
    lex.order = TRUE
  )
  pieces <- lapply(split(seq_len(nrow(data)), group_id), function(index) {
    group_values <- data[index[1L], grouping, drop = FALSE]
    values <- lapply(metric_names, function(metric) {
      x <- data[[metric]][index]
      c(
        mean = mean(x),
        sd = stats::sd(x),
        se = stats::sd(x) / sqrt(length(x))
      )
    })
    names(values) <- metric_names
    wide <- unlist(values)
    result <- cbind(
      data.frame(analysis_protocol = analysis_protocol),
      group_values,
      data.frame(
        runs = length(index),
        as.list(wide),
        check.names = FALSE
      )
    )
    names(result) <- gsub(".", "_", names(result), fixed = TRUE)
    result
  })
  result <- do.call(rbind, pieces)
  rownames(result) <- NULL
  result
}

summary_all <- summarize_metrics(
  metrics,
  c("method", "nominal_fdr")
)
summary_all$scope <- "all scenarios"
summary_baseline <- summarize_metrics(
  metrics[metrics$is_baseline, , drop = FALSE],
  c("method", "nominal_fdr")
)
summary_baseline$scope <- "four-replicate baseline"
summary_table <- rbind(summary_all, summary_baseline)
write.csv(
  summary_table,
  file.path(
    "data", "derived", "crispulator_facs_f1_by_fdr_summary.csv"
  ),
  row.names = FALSE
)

replicate_summary <- summarize_metrics(
  metrics,
  c("replicates", "method", "nominal_fdr")
)
write.csv(
  replicate_summary,
  file.path(
    "data", "derived",
    "crispulator_facs_f1_by_fdr_replicate_summary.csv"
  ),
  row.names = FALSE
)

method_colours <- c(
  `BARCS-original` = "#0072B2",
  `BARCS-NORM` = "#7A3E9D",
  `BARCS-partial` = "#009E73",
  `BARCS-EB` = "#D55E00",
  `MAGeCK-MLE` = "#6A3D9A",
  `edgeR-QL` = "#CC79A7",
  DESeq2 = "#56B4E9",
  `limma-voom` = "#666666"
)
method_pch <- c(16, 17, 15, 18, 19, 8, 4, 3)
x_positions <- -log10(thresholds)
x_labels <- format(thresholds, trim = TRUE, scientific = FALSE)

draw_panel <- function(scope, metric, title, ylim) {
  panel <- summary_table[
    summary_table$scope == scope,
    ,
    drop = FALSE
  ]
  plot(
    NA,
    xlim = range(x_positions),
    ylim = ylim,
    xaxt = "n",
    xlab = "Nominal gene FDR",
    ylab = if (metric == "f1") "F1 score" else "Realized FDP",
    main = title,
    bty = "l"
  )
  axis(1, at = x_positions, labels = x_labels)
  if (metric == "realized_fdp") {
    lines(
      x_positions,
      thresholds,
      lty = 2,
      lwd = 1.5,
      col = "#333333"
    )
  }
  for (method_index in seq_along(methods)) {
    method <- methods[method_index]
    method_panel <- panel[panel$method == method, ]
    method_panel <- method_panel[
      match(thresholds, method_panel$nominal_fdr),
      ,
      drop = FALSE
    ]
    lines(
      x_positions,
      method_panel[[paste0(metric, "_mean")]],
      type = "o",
      pch = method_pch[method_index],
      lwd = 2,
      col = method_colours[[method]]
    )
  }
}

pdf(
  file.path("figures", "crispulator_facs_f1_by_fdr.pdf"),
  width = 10.5,
  height = 8,
  useDingbats = FALSE
)
layout(
  matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE),
  heights = c(1, 1, 0.22)
)
par(mar = c(4.5, 4.3, 2.8, 1))
draw_panel("all scenarios", "f1", "A  F1: all 50 runs", c(0, 1))
draw_panel(
  "four-replicate baseline",
  "f1",
  "B  F1: baseline five seeds",
  c(0, 1)
)
draw_panel(
  "all scenarios",
  "realized_fdp",
  "C  Realized FDP: all 50 runs",
  c(0, 0.60)
)
draw_panel(
  "four-replicate baseline",
  "realized_fdp",
  "D  Realized FDP: baseline five seeds",
  c(0, 0.60)
)
par(mar = c(0, 0, 0, 0))
plot.new()
legend(
  "center",
  legend = methods,
  col = unname(method_colours[methods]),
  pch = method_pch,
  lty = 1,
  lwd = 2,
  horiz = TRUE,
  bty = "n",
  cex = 0.78
)
dev.off()

print(summary_table)
