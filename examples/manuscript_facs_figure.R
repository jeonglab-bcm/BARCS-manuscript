#!/usr/bin/env Rscript

# Main-text figure for the IL2RA ordered-bin benchmark.
#
# The values are read from the reproducible benchmark outputs. Published
# Waterbear and MAUDE aggregates are kept separate from independently rerun
# methods and are used only in panels where the deposited paper reports them.

source(file.path("R", "method_palette.R"))

metrics_path <- file.path(
  "data", "derived", "waterbear_facs_external_head_to_head_metrics.csv"
)
null_path <- file.path(
  "data", "derived", "waterbear_facs_three_method_null_calibration.csv"
)
output_path <- file.path("figures", "waterbear_facs_main.pdf")

metrics <- read.csv(metrics_path, check.names = FALSE)
null_metrics <- read.csv(null_path, check.names = FALSE)

method_ids <- c("BARCS-original", "MAGeCK-MLE", "Waterbear", "MAUDE")
method_labels <- c("BARCS", "MAGeCK-MLE", "Waterbear", "MAUDE")
method_rows <- match(method_ids, metrics$method)
if (anyNA(method_rows)) {
  stop("The FACS benchmark output is missing a main-text method.")
}
metrics <- metrics[method_rows, ]

method_colours <- c(
  barcs_method_colours[["BARCS"]],
  barcs_method_colours[["MAGeCK"]],
  barcs_method_colours[["Waterbear"]],
  barcs_method_colours[["MAUDE"]]
)

pdf(output_path, width = 10.5, height = 3.7, useDingbats = FALSE)
layout(matrix(1:3, nrow = 1), widths = c(1, 1, 0.95))
par(family = "sans", las = 1, cex.axis = 0.88, cex.lab = 0.95)

par(mar = c(5.2, 4.1, 2.7, 0.7))
recovery <- metrics$directionally_validated_recovered
recovery_positions <- barplot(
  recovery,
  names.arg = method_labels,
  col = method_colours,
  border = NA,
  ylim = c(0, 27),
  ylab = "Validated regulators recovered",
  main = "A  Directional recovery",
  las = 2,
  cex.names = 0.78
)
abline(h = 26, lty = 3, col = "#777777")
text(
  recovery_positions,
  recovery,
  labels = paste0(recovery, "/26"),
  pos = 3,
  cex = 0.78
)

par(mar = c(5.2, 4.1, 2.7, 0.7))
discoveries <- metrics$discoveries_at_fdr_0_10
call_positions <- barplot(
  discoveries,
  names.arg = method_labels,
  col = method_colours,
  border = NA,
  log = "y",
  ylim = c(10, 600),
  yaxt = "n",
  ylab = "Screen discoveries (log scale)",
  main = "B  Calls at method threshold",
  las = 2,
  cex.names = 0.78
)
axis(2, at = c(10, 30, 100, 300, 600), labels = c(10, 30, 100, 300, 600))
text(
  call_positions,
  discoveries,
  labels = discoveries,
  pos = 3,
  cex = 0.78
)

par(mar = c(5.2, 4.3, 2.7, 0.7))
null_rate <- null_metrics$fraction_p_below_0_05[
  match(c("raw", "control_calibrated"), null_metrics$guide_statistic)
]
null_positions <- barplot(
  null_rate,
  names.arg = c("Raw", "Control-\ncalibrated"),
  col = c("#9ECAE1", barcs_method_colours[["BARCS"]]),
  border = NA,
  ylim = c(0, 0.15),
  ylab = "Non-targeting guides with p < 0.05",
  main = "C  Empirical null",
  las = 1,
  cex.names = 0.82
)
abline(h = 0.05, lty = 2, lwd = 1.4, col = "#555555")
text(
  null_positions,
  null_rate,
  labels = sprintf("%.3f", null_rate),
  pos = 3,
  cex = 0.8
)
text(
  par("usr")[2], 0.052, "nominal 0.05",
  adj = c(1, 0), cex = 0.72, col = "#555555"
)

dev.off()
