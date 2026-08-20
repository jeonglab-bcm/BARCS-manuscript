#!/usr/bin/env Rscript

# Figure for the simCRISPR interaction benchmark.
#
# Split from the benchmark script so the plot can be redrawn without refitting.
# It reads only the committed tables in `data/derived`, so it runs in seconds:
#
#     Rscript examples/simcrispr_interaction_figure.R
#
# The experiment varies one thing: what BARCS divides guide counts by. Three of
# the four series are therefore the same estimator run three times. MAGeCK-MLE
# has no equivalent choice to sweep and appears once, as a reference. An
# earlier draft put all four in one legend headed "denominator", which mixed
# two different kinds of thing; every series is now labelled with its method
# first so the grouping is unambiguous.

options(stringsAsFactors = FALSE)
source(file.path("R", "method_palette.R"))
derived <- file.path("data", "derived")
scan_path <- file.path(derived, "simcrispr_interaction_scan_summary.csv")
null_path <- file.path(derived, "simcrispr_interaction_null_estimates.csv")
for (path in c(scan_path, null_path)) {
  if (!file.exists(path)) {
    stop(
      "Missing ", path,
      ". Run examples/simcrispr_interaction_benchmark.R first.",
      call. = FALSE
    )
  }
}
scan_summary <- read.csv(scan_path, check.names = FALSE)
null_estimates <- read.csv(null_path, check.names = FALSE)

# Series definitions. The label carries the method name because that is the
# distinction readers were losing.
#
# Colours come from R/method_palette.R, which gives one base colour per method
# and separates variants of the same method by point shape, line type, or fill.
# The three BARCS denominators are one estimator run three times, so they stay
# within the BARCS hue. The previous hard-coded vector spent #D55E00 -- MAGeCK's
# base colour -- on a BARCS series and drew MAGeCK-MLE in #6A3D9A, so the same
# method was orange in Figure 2 and violet here.
barcs_colour <- barcs_method_colours[["BARCS-moderated"]]
mageck_colour <- barcs_method_colours[["MAGeCK-MLE"]]
# The non-targeting and safe-harbor curves in panel B very nearly coincide, so
# drawing all three denominators in one flat colour makes them inseparable even
# with different line types. They instead take three tints of the BARCS blue:
# one hue, so the reader still sees one method, with enough lightness contrast
# to follow overlapping curves.
barcs_tints <- c(barcs_colour, "#5FA8D3", "#00446E")
series <- list(
  list(
    key = "BARCS-moderated [library]",
    label = "BARCS, divided by whole library",
    colour = barcs_tints[1L], pch = 16, lty = 1
  ),
  list(
    key = "BARCS-moderated [non-targeting]",
    label = "BARCS, divided by non-targeting controls",
    colour = barcs_tints[2L], pch = 17, lty = 2
  ),
  list(
    key = "BARCS-moderated [safe-harbor]",
    label = "BARCS, divided by safe-harbor controls",
    colour = barcs_tints[3L], pch = 15, lty = 3
  ),
  list(
    key = "MAGeCK-MLE [non-targeting]",
    label = "MAGeCK-MLE (reference; no such choice; B only)",
    colour = mageck_colour, pch = 8, lty = 1
  )
)
barcs_series <- series[1:3]
thresholds <- sort(unique(scan_summary$nominal_fdr))

pdf(
  file.path("figures", "simcrispr_interaction.pdf"),
  width = 10, height = 4.8, useDingbats = FALSE
)
layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(1, 0.26))

# ---- Panel A: the bias ------------------------------------------------------
# Guides whose true interaction is exactly zero. Anything away from the dashed
# line is the denominator inventing an effect that is not there.
par(mar = c(4.6, 12.5, 3.4, 1), bty = "l")
values <- lapply(barcs_series, function(entry) {
  null_estimates$estimate[null_estimates$method == entry$key]
})
names(values) <- vapply(barcs_series, function(entry) {
  sub("^BARCS, divided by ", "", entry$label)
}, character(1))
# The fixed c(-0.32, 0.32) cut the whole-library upper whisker, which reaches
# 0.370 -- the panel exists to show that this denominator invents an effect, and
# the axis was hiding how large it gets. Limits now come from the whiskers, and
# the count of suppressed outliers is stated rather than left implicit.
whisker_range <- range(vapply(
  values, function(v) boxplot.stats(v)$stats[c(1L, 5L)], numeric(2)
))
value_limit <- c(-1, 1) * (1.06 * max(abs(whisker_range)))
suppressed <- sum(vapply(
  values, function(v) length(boxplot.stats(v)$out), integer(1)
))
boxplot(
  rev(values), horizontal = TRUE, las = 1, outline = FALSE,
  col = rev(vapply(barcs_series, function(entry) entry$colour, character(1))),
  border = "#333333", ylim = value_limit,
  xlab = "Estimated interaction (truth is exactly 0)",
  main = "A  Guides that did nothing", cex.main = 1.1
)
abline(v = 0, lty = 2, lwd = 1.6, col = "#333333")
mtext(
  sprintf(
    "whiskers to 1.5 IQR; %d outlier%s not drawn",
    suppressed, if (suppressed == 1L) "" else "s"
  ),
  side = 3, line = 0.15, adj = 1, cex = 0.62, col = "#666666"
)

# ---- Panel B: what the bias costs ------------------------------------------
par(mar = c(4.6, 4.6, 3.4, 1), pty = "s", bty = "l")
x_at <- -log10(thresholds)
# 0.2, 0.1 and 0.05 sit 0.30 and 0.30 apart on the -log10 axis, so labelling all
# six thresholds ran "0.2", "0.1" and "0.05" into each other. Every threshold
# keeps a minor tick; only a non-colliding subset is labelled, matching the
# treatment in examples/manuscript_crispulator_figure.R.
labelled_thresholds <- c(0.2, 0.01, 1e-3, 1e-4)
x_text <- as.expression(lapply(labelled_thresholds, function(value) {
  exponent <- round(log10(value))
  if (isTRUE(all.equal(value, 10^exponent)) && exponent <= -3) {
    bquote(10^.(exponent))
  } else {
    format(value, trim = TRUE, scientific = FALSE)
  }
}))
plot(
  NA, xlim = range(x_at), ylim = c(0, 1), xaxt = "n",
  xlab = "Nominal gene FDR", ylab = "F1 score",
  main = "B  Real interactions recovered", bty = "l", cex.main = 1.1
)
axis(1, at = x_at, labels = FALSE, tcl = -0.2)
axis(1, at = -log10(labelled_thresholds), labels = x_text)
for (entry in series) {
  rows <- scan_summary[
    scan_summary$method == entry$key & scan_summary$metric == "f1", ,
    drop = FALSE
  ]
  rows <- rows[match(thresholds, rows$nominal_fdr), , drop = FALSE]
  if (anyNA(rows$mean)) {
    stop(
      entry$key, " is missing ", sum(is.na(rows$mean)), " of ",
      length(thresholds),
      " thresholds; the curve would be drawn with unmarked gaps.",
      call. = FALSE
    )
  }
  lines(x_at, rows$mean, type = "o", pch = entry$pch, lty = entry$lty,
        lwd = 2, col = entry$colour)
}

par(mar = c(0, 0, 0, 0), pty = "m")
plot.new()
legend(
  "center",
  legend = vapply(series, function(entry) entry$label, character(1)),
  col = vapply(series, function(entry) entry$colour, character(1)),
  pch = vapply(series, function(entry) entry$pch, numeric(1)),
  lty = vapply(series, function(entry) entry$lty, numeric(1)),
  lwd = 2, ncol = 2, bty = "n", cex = 1.15
)
dev.off()

cat("Wrote figures/simcrispr_interaction.pdf\n")
