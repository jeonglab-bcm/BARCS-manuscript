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
if (!file.exists(summary_path)) {
  stop(
    "Missing ", summary_path,
    ". Run examples/crispulator_facs_benchmark.R first.",
    call. = FALSE
  )
}
summary_table <- read.csv(summary_path, check.names = FALSE)
design <- "low-bulk-high"
if (!design %in% summary_table$design) {
  stop(
    "Design ", design, " is absent from ", summary_path,
    ". Available: ", paste(sort(unique(summary_table$design)), collapse = ", "),
    call. = FALSE
  )
}
summary_table <- summary_table[summary_table$design == design, ]

methods <- c(
  "BARCS-original", "BARCS-moderated", "MAGeCK-MLE", "CRISPhieRmix"
)
method_labels <- c(
  "BARCS-ST", "BARCS-MOD", "MAGeCK-MLE", "CRISPhieRmix"
)
# R/method_palette.R is the single source of truth: one base colour per method,
# with variants of the same method sharing it and separating by point shape and
# line type. BARCS-ST and BARCS-MOD are two settings of one estimator, so they
# share the BARCS blue. The previous hard-coded vector gave BARCS-ST #56B4E9,
# which is CB2's base colour, and sourced the palette without using it.
colours <- barcs_method_colours[
  c("BARCS-ST", "BARCS-MOD", "MAGeCK-MLE", "CRISPhieRmix")
]
if (anyNA(colours)) {
  stop(
    "R/method_palette.R does not define: ",
    paste(
      c("BARCS-ST", "BARCS-MOD", "MAGeCK-MLE", "CRISPhieRmix")[is.na(colours)],
      collapse = ", "
    ),
    call. = FALSE
  )
}
symbols <- c(16, 17, 15, 8)
line_types <- c(1, 2, 1, 1)
names(colours) <- names(symbols) <- names(line_types) <- methods

thresholds <- sort(unique(summary_table$nominal_fdr))
x <- -log10(thresholds)
# 0.20 and 0.10 are 0.30 apart on the -log10 axis while every other labelled
# pair is 1.00 apart, so labelling both collides. All 13 thresholds still get a
# minor tick.
labelled_thresholds <- c(0.20, 0.01, 1e-3, 1e-4, 1e-5, 1e-6)
axis_labels <- as.expression(lapply(labelled_thresholds, function(value) {
  exponent <- round(log10(value))
  if (isTRUE(all.equal(value, 10^exponent)) && exponent <= -3) {
    bquote(10^.(exponent))
  } else {
    format(value, trim = TRUE, scientific = FALSE)
  }
}))

# Panels B and D previously used a fixed c(0, 0.30). CRISPhieRmix reaches a
# realized FDP of 0.358 at MOI 0.20 and 0.379 at MOI 0.30, so its worst point
# fell outside the panel and was drawn off the top edge with nothing to mark it
# -- the clipping hid a competitor's largest error. Both metrics now take limits
# from the data, and `draw_panel` refuses to clip silently.
metric_limits <- list(
  f1 = c(0, max(0.95, 1.02 * max(summary_table$f1))),
  realized_fdp = c(0, max(0.30, 1.05 * max(summary_table$realized_fdp)))
)

draw_panel <- function(moi, metric, title) {
  panel <- summary_table[summary_table$moi == moi, ]
  y_limit <- metric_limits[[metric]]
  if (max(panel[[metric]]) > y_limit[2L] || min(panel[[metric]]) < y_limit[1L]) {
    stop(
      "Panel '", title, "' would clip ", metric, ": data range [",
      paste(signif(range(panel[[metric]]), 4), collapse = ", "),
      "] against limits [", paste(y_limit, collapse = ", "), "]",
      call. = FALSE
    )
  }
  plot(
    NA,
    xlim = range(x),
    ylim = y_limit,
    xaxt = "n",
    xlab = "Nominal gene FDR",
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
    if (anyNA(rows[[metric]])) {
      stop(
        method, " is missing ", sum(is.na(rows[[metric]])), " of ",
        length(thresholds), " thresholds at MOI ", moi,
        "; the curve would be drawn with unmarked gaps.",
        call. = FALSE
      )
    }
    lines(
      x,
      rows[[metric]],
      type = "o",
      pch = symbols[[method]],
      lty = line_types[[method]],
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
par(mar = c(4.4, 4.2, 2.5, 0.8), family = "sans", pty = "s")
draw_panel(0.20, "f1", "A  F1: MOI 0.20")
draw_panel(0.20, "realized_fdp", "B  Realized FDP: MOI 0.20")
draw_panel(0.30, "f1", "C  F1: MOI 0.30")
draw_panel(0.30, "realized_fdp", "D  Realized FDP: MOI 0.30")
par(mar = c(0, 0, 0, 0), pty = "m")
plot.new()
# The dashed diagonal in B and D was drawn but never named, so the one line a
# reader must understand to judge FDP control was unlabelled.
legend(
  "center",
  legend = c(method_labels, "Realized FDP = nominal (B, D)"),
  col = c(unname(colours[methods]), "#444444"),
  pch = c(unname(symbols[methods]), NA),
  lty = c(unname(line_types[methods]), 2),
  lwd = 2,
  horiz = TRUE,
  bty = "n",
  cex = 1.1
)
dev.off()
