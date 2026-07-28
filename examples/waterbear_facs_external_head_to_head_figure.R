#!/usr/bin/env Rscript

# Figure for the GSE242880 IL2RA head-to-head comparison.
#
# Split from the benchmark script so the plot can be redrawn without refitting.
# It reads only the committed table in `data/derived`, so it runs in seconds:
#
#     Rscript examples/waterbear_facs_external_head_to_head_figure.R
#
# BARCS-NORM is deliberately not drawn. It makes zero calls at FDR 0.10 on
# this screen, so its bars in panels A and B are empty and its F1 and balanced
# accuracy in panel C are 0 and 0.5 by construction. Three of the eight columns
# therefore carried no information while costing the reader a category to
# track. The method is retained in the metrics CSV written by the benchmark
# script and discussed in `docs/barcs-external-method-comparison.md`; it is
# the figure, not the record, that drops it.

options(stringsAsFactors = FALSE)
derived <- file.path("data", "derived")
metrics_path <- file.path(
  derived, "waterbear_facs_external_head_to_head_metrics.csv"
)
if (!file.exists(metrics_path)) {
  stop(
    "Missing ", metrics_path,
    ". Run examples/waterbear_facs_external_head_to_head.R first.",
    call. = FALSE
  )
}
metrics <- read.csv(metrics_path, check.names = FALSE)

omitted <- "BARCS-NORM"
metrics <- metrics[!(metrics$method %in% omitted), , drop = FALSE]

# Display order: the BARCS gene statistics, then the two independently rerun
# MAGeCK analyses, then the two literature aggregates that have no per-gene
# output and so cannot enter panel C.
method_order <- c(
  "BARCS-original", "BARCS-partial", "BARCS-EB",
  "MAGeCK-MLE", "MAGeCK-test", "Waterbear", "MAUDE"
)
if (!setequal(metrics$method, method_order)) {
  stop(
    "Unexpected method set in ", metrics_path, ": ",
    paste(sort(setdiff(metrics$method, method_order)), collapse = ", "),
    call. = FALSE
  )
}
metrics <- metrics[match(method_order, metrics$method), , drop = FALSE]
rerun_metrics <- metrics[
  metrics$result_type == "independent rerun", ,
  drop = FALSE
]

# Colours are inherited unchanged from the two CRISPulator gene-statistic
# figures so that one BARCS variant keeps one colour across the family.
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

# ---- Panel A: how much validated biology each method returns ---------------
par(mar = c(7.6, 4.2, 2.7, 0.8))
positions <- barplot(
  metrics$directionally_validated_recovered,
  names.arg = metrics$method,
  col = unname(method_colours[metrics$method]),
  border = NA,
  las = 2,
  cex.names = 0.72,
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
  cex = 0.76
)

# ---- Panel B: what it cost in calls to get there ---------------------------
par(mar = c(7.6, 4.2, 2.7, 0.8))
positions <- barplot(
  metrics$discoveries_at_fdr_0_10,
  names.arg = metrics$method,
  col = unname(method_colours[metrics$method]),
  border = NA,
  las = 2,
  cex.names = 0.72,
  ylim = c(0, 440),
  ylab = "Screen discoveries",
  main = "B  Calls at reported threshold"
)
text(
  positions,
  metrics$discoveries_at_fdr_0_10,
  labels = metrics$discoveries_at_fdr_0_10,
  pos = 3,
  cex = 0.76
)

# ---- Panel C: recovery and precision read together -------------------------
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
  cex.names = 0.72,
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
  cex = 0.68,
  horiz = TRUE,
  inset = 0.01
)
dev.off()

cat(
  "Wrote figures/waterbear_facs_external_head_to_head.pdf without ",
  paste(omitted, collapse = ", "), "\n",
  sep = ""
)
