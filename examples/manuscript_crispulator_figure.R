#!/usr/bin/env Rscript

# Redraw the genome-scale CRISPulator result from its derived summary table.
# This script has no simulator or method-package dependency.

source(file.path("R", "method_palette.R"))

summary_path <- file.path(
  "data", "derived", "crispulator_facs_moi_10k_f1_by_fdr_summary.csv"
)
output_path <- file.path(
  "figures", "crispulator_facs_moi_10k_f1_by_fdr.pdf"
)
summary_table <- read.csv(summary_path, check.names = FALSE)
summary_table <- summary_table[summary_table$design == "low-bulk-high", ]

methods <- c(
  "BARCS-original", "BARCS-moderated", "MAGeCK-MLE", "CRISPhieRmix"
)
method_labels <- c(
  "BARCS-ST", "BARCS-MOD", "MAGeCK-MLE", "CRISPhieRmix"
)
colours <- c("#56B4E9", "#0072B2", "#D55E00", "#7A3E9D")
symbols <- c(16, 17, 15, 8)
names(colours) <- names(symbols) <- methods

thresholds <- sort(unique(summary_table$nominal_fdr))
x <- -log10(thresholds)
labelled_thresholds <- c(0.20, 0.10, 0.01, 1e-3, 1e-4, 1e-5, 1e-6)
axis_labels <- as.expression(lapply(labelled_thresholds, function(value) {
  exponent <- round(log10(value))
  if (isTRUE(all.equal(value, 10^exponent)) && exponent <= -3) {
    bquote(10^.(exponent))
  } else {
    format(value, trim = TRUE, scientific = FALSE)
  }
}))

draw_panel <- function(moi, metric, title, y_limit) {
  panel <- summary_table[summary_table$moi == moi, ]
  plot(
    NA,
    xlim = range(x),
    ylim = y_limit,
    xaxt = "n",
    xlab = "Requested gene FDR",
    ylab = if (metric == "f1") "F1" else "Realized FDP",
    main = title,
    bty = "l"
  )
  axis(1, at = x, labels = FALSE, tcl = -0.2)
  axis(1, at = -log10(labelled_thresholds), labels = axis_labels)
  if (metric == "realized_fdp") {
    lines(x, thresholds, lty = 2, lwd = 1.4, col = "#444444")
  }
  for (method in methods) {
    rows <- panel[panel$method == method, ]
    rows <- rows[match(thresholds, rows$nominal_fdr), ]
    lines(
      x,
      rows[[metric]],
      type = "o",
      pch = symbols[[method]],
      lwd = 2,
      col = colours[[method]]
    )
  }
}

pdf(output_path, width = 10.5, height = 7.5, useDingbats = FALSE)
layout(
  matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE),
  heights = c(1, 1, 0.16)
)
par(mar = c(4.4, 4.2, 2.5, 0.8), family = "sans")
draw_panel(0.20, "f1", "A  F1: MOI 0.20", c(0, 0.95))
draw_panel(0.20, "realized_fdp", "B  Realized FDP: MOI 0.20", c(0, 0.30))
draw_panel(0.30, "f1", "C  F1: MOI 0.30", c(0, 0.95))
draw_panel(0.30, "realized_fdp", "D  Realized FDP: MOI 0.30", c(0, 0.30))
par(mar = c(0, 0, 0, 0))
plot.new()
legend(
  "center",
  legend = method_labels,
  col = unname(colours[methods]),
  pch = unname(symbols[methods]),
  lty = 1,
  lwd = 2,
  horiz = TRUE,
  bty = "n",
  cex = 0.86
)
dev.off()
