#!/usr/bin/env Rscript

# Aggregation stage of the CRISPulator head-to-head comparison.
#
# `examples/crispulator_facs_external_head_to_head.R` builds the per-run
# metric table from the gene-result files under the ignored `results/` tree
# and then sources this file. Because the per-run table is committed under
# `data/derived/`, this stage can also be re-run on its own to rebuild the
# summary tables and the figure without the large result directories:
#
#     Rscript examples/crispulator_facs_external_head_to_head_aggregate.R
#
# Analysis scope. The headline comparison is restricted to the settings with
# at least `min_replicates` independent screen replicates. The one-replicate
# setting is a prespecified diagnostic boundary rather than a supported
# operating point: at R = 1 the low--bulk--high design leaves a single
# residual degree of freedom per guide, so the beta-binomial guide test has
# no usable sample-replicate reference distribution and BARCS-original makes
# no calls at all. Pooling that setting into the headline mean reports a
# boundary failure as if it were average performance. It is still summarised
# here, separately and under its own scope label, because the comparison
# against the count models at R = 1 is itself a result.

options(stringsAsFactors = FALSE)

if (!exists("analysis_protocol", inherits = FALSE)) {
  analysis_protocol <- "barcs-external-headtohead-v1"
}
if (!exists("expected_methods", inherits = FALSE)) {
  expected_methods <- c(
    "BARCS-original", "BARCS-NORM", "BARCS-partial", "BARCS-EB",
    "MAGeCK-MLE", "edgeR-QL", "DESeq2", "limma-voom"
  )
}

metrics_path <- file.path(
  "data", "derived",
  "crispulator_facs_external_head_to_head_metrics.csv"
)
if (!exists("metrics", inherits = FALSE)) {
  if (!file.exists(metrics_path)) {
    stop(
      "Per-run metrics are unavailable. Run ",
      "examples/crispulator_facs_external_head_to_head.R first."
    )
  }
  metrics <- read.csv(metrics_path)
}
stopifnot(
  identical(sort(unique(metrics$method)), sort(expected_methods)),
  all(metrics$analysis_protocol == analysis_protocol),
  "replicates" %in% names(metrics)
)

min_replicates <- 2L
primary_scope <- "multi-replicate settings"
boundary_scope <- "one-replicate boundary"
scope_of <- function(replicates) {
  ifelse(replicates >= min_replicates, primary_scope, boundary_scope)
}
metrics$scope <- scope_of(metrics$replicates)
scopes <- c(primary_scope, boundary_scope)
scopes <- scopes[scopes %in% metrics$scope]
stopifnot(primary_scope %in% scopes)

summary_metrics <- c(
  "auroc", "average_precision", "directional_recall_fdr_0_10",
  "empirical_fdp_fdr_0_10", "f1_fdr_0_10",
  "negative_control_p_below_0_05"
)

summary_rows <- list()
summary_index <- 0L
for (scope in scopes) {
  scope_data <- metrics[metrics$scope == scope, , drop = FALSE]
  for (method in expected_methods) {
    method_data <- scope_data[scope_data$method == method, , drop = FALSE]
    for (metric in summary_metrics) {
      values <- method_data[[metric]]
      summary_index <- summary_index + 1L
      summary_rows[[summary_index]] <- data.frame(
        analysis_protocol = analysis_protocol,
        scope = scope,
        method = method,
        metric = metric,
        runs = length(values),
        mean = mean(values),
        sd = sd(values),
        se = sd(values) / sqrt(length(values)),
        minimum = min(values),
        maximum = max(values)
      )
    }
  }
}
summary_table <- do.call(rbind, summary_rows)
write.csv(
  summary_table,
  file.path(
    "data", "derived",
    "crispulator_facs_external_head_to_head_summary.csv"
  ),
  row.names = FALSE
)

paired_rows <- list()
paired_index <- 0L
references <- c("BARCS-original", "BARCS-EB")
for (scope in scopes) {
  scope_data <- metrics[metrics$scope == scope, , drop = FALSE]
  for (reference in references) {
    reference_data <- scope_data[
      scope_data$method == reference,
      c("scenario_id", "seed", summary_metrics)
    ]
    for (method in setdiff(expected_methods, reference)) {
      method_data <- scope_data[
        scope_data$method == method,
        c("scenario_id", "seed", summary_metrics)
      ]
      paired <- merge(
        method_data,
        reference_data,
        by = c("scenario_id", "seed"),
        suffixes = c("_method", "_reference")
      )
      for (metric in summary_metrics) {
        difference <- paired[[paste0(metric, "_method")]] -
          paired[[paste0(metric, "_reference")]]
        half_width <- qt(0.975, df = length(difference) - 1L) *
          sd(difference) / sqrt(length(difference))
        paired_index <- paired_index + 1L
        paired_rows[[paired_index]] <- data.frame(
          analysis_protocol = analysis_protocol,
          scope = scope,
          method = method,
          reference = reference,
          metric = metric,
          runs = length(difference),
          mean_difference = mean(difference),
          sd_difference = sd(difference),
          lower_95 = mean(difference) - half_width,
          upper_95 = mean(difference) + half_width
        )
      }
    }
  }
}
paired_table <- do.call(rbind, paired_rows)
write.csv(
  paired_table,
  file.path(
    "data", "derived",
    "crispulator_facs_external_head_to_head_paired.csv"
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
figure_metrics <- c(
  average_precision = "Average precision",
  directional_recall_fdr_0_10 = "Directional recall",
  empirical_fdp_fdr_0_10 = "Realized FDP"
)
pdf(
  file.path("figures", "crispulator_facs_external_head_to_head.pdf"),
  width = 11,
  height = 4.4,
  useDingbats = FALSE
)
layout(matrix(1:3, nrow = 1))
for (panel_index in seq_along(figure_metrics)) {
  metric <- names(figure_metrics)[panel_index]
  panel <- summary_table[
    summary_table$scope == primary_scope &
      summary_table$metric == metric,
    ,
    drop = FALSE
  ]
  panel <- panel[match(expected_methods, panel$method), ]
  upper <- panel$mean + panel$se
  ylim <- c(0, max(upper) * 1.12)
  if (metric != "empirical_fdp_fdr_0_10") {
    ylim <- c(0, 1)
  }
  best_index <- if (metric == "empirical_fdp_fdr_0_10") {
    which.min(panel$mean)
  } else {
    which.max(panel$mean)
  }
  borders <- rep(NA_character_, nrow(panel))
  borders[best_index] <- "#111111"
  par(mar = c(8.2, 4.2, 2.6, 0.8))
  positions <- barplot(
    panel$mean,
    names.arg = panel$method,
    col = unname(method_colours[panel$method]),
    border = borders,
    lwd = 2,
    las = 2,
    cex.names = 0.67,
    ylim = ylim,
    ylab = figure_metrics[[metric]],
    main = paste0(
      LETTERS[panel_index], "  ",
      figure_metrics[[metric]],
      if (metric == "empirical_fdp_fdr_0_10") {
        " (lower is safer)"
      } else {
        ""
      }
    )
  )
  arrows(
    positions,
    panel$mean - panel$se,
    positions,
    panel$mean + panel$se,
    angle = 90,
    code = 3,
    length = 0.035
  )
  if (metric == "empirical_fdp_fdr_0_10") {
    abline(h = 0.10, lty = 2, col = "#555555")
  }
}
dev.off()

print(summary_table)
print(paired_table)
