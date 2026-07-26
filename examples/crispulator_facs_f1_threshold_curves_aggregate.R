#!/usr/bin/env Rscript

# Aggregation stage of the CRISPulator F1-versus-threshold curves.
#
# `examples/crispulator_facs_f1_threshold_curves.R` builds the per-run
# threshold table from the gene-result files under the ignored `results/`
# tree and then sources this file. Because the per-run table is committed
# under `data/derived/`, this stage can also be re-run on its own:
#
#     Rscript examples/crispulator_facs_f1_threshold_curves_aggregate.R
#
# Analysis scope. The pooled curve covers the multi-replicate settings only.
# At R = 1 the low--bulk--high design leaves one residual degree of freedom
# per guide, BARCS-original makes no calls at any threshold, and the count
# models run at realized FDP far above nominal, so pooling that setting into
# the headline curve mixes a boundary failure with ordinary operating points.
# It is drawn separately instead.

options(stringsAsFactors = FALSE)

if (!exists("analysis_protocol", inherits = FALSE)) {
  analysis_protocol <- "crispulator-f1-threshold-curves-v1"
}
if (!exists("thresholds", inherits = FALSE)) {
  thresholds <- c(0.10, 0.05, 0.01, 0.005, 0.001)
}
if (!exists("methods", inherits = FALSE)) {
  methods <- c(
    "BARCS-original", "BARCS-NORM", "BARCS-partial", "BARCS-EB",
    "MAGeCK-MLE", "edgeR-QL", "DESeq2", "limma-voom"
  )
}

metrics_path <- file.path(
  "data", "derived", "crispulator_facs_f1_by_fdr.csv"
)
if (!exists("metrics", inherits = FALSE)) {
  if (!file.exists(metrics_path)) {
    stop(
      "Per-run threshold metrics are unavailable. Run ",
      "examples/crispulator_facs_f1_threshold_curves.R first."
    )
  }
  metrics <- read.csv(metrics_path)
}
stopifnot(
  setequal(unique(metrics$method), methods),
  setequal(unique(metrics$nominal_fdr), thresholds),
  "replicates" %in% names(metrics)
)

min_replicates <- 2L
primary_scope <- "multi-replicate settings"
boundary_scope <- "one-replicate boundary"

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

multi_replicate <- metrics[
  metrics$replicates >= min_replicates, ,
  drop = FALSE
]
one_replicate <- metrics[
  metrics$replicates < min_replicates, ,
  drop = FALSE
]
stopifnot(nrow(multi_replicate) > 0L)

summary_primary <- summarize_metrics(
  multi_replicate,
  c("method", "nominal_fdr")
)
summary_primary$scope <- primary_scope
summary_baseline <- summarize_metrics(
  metrics[metrics$is_baseline, , drop = FALSE],
  c("method", "nominal_fdr")
)
summary_baseline$scope <- "four-replicate baseline"
summary_table <- rbind(summary_primary, summary_baseline)
if (nrow(one_replicate) > 0L) {
  summary_boundary <- summarize_metrics(
    one_replicate,
    c("method", "nominal_fdr")
  )
  summary_boundary$scope <- boundary_scope
  summary_table <- rbind(summary_table, summary_boundary)
}
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
draw_panel(
  primary_scope,
  "f1",
  "A  F1: 45 multi-replicate runs",
  c(0, 1)
)
draw_panel(
  "four-replicate baseline",
  "f1",
  "B  F1: baseline five seeds",
  c(0, 1)
)
draw_panel(
  primary_scope,
  "realized_fdp",
  "C  Realized FDP: 45 multi-replicate runs",
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
