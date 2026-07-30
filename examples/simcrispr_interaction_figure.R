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
series <- list(
  list(
    key = "BARCS-moderated [library]",
    label = "BARCS, divided by whole library",
    colour = "#0072B2", pch = 16
  ),
  list(
    key = "BARCS-moderated [non-targeting]",
    label = "BARCS, divided by non-targeting controls",
    colour = "#D55E00", pch = 17
  ),
  list(
    key = "BARCS-moderated [safe-harbor]",
    label = "BARCS, divided by safe-harbor controls",
    colour = "#009E73", pch = 15
  ),
  list(
    key = "MAGeCK-MLE [non-targeting]",
    label = "MAGeCK-MLE (reference; no such choice)",
    colour = "#6A3D9A", pch = 8
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
par(mar = c(4.6, 12.5, 3.0, 1))
values <- lapply(barcs_series, function(entry) {
  null_estimates$estimate[null_estimates$method == entry$key]
})
names(values) <- vapply(barcs_series, function(entry) {
  sub("^BARCS, divided by ", "", entry$label)
}, character(1))
boxplot(
  rev(values), horizontal = TRUE, las = 1, outline = FALSE,
  col = rev(vapply(barcs_series, function(entry) entry$colour, character(1))),
  border = "#333333", ylim = c(-0.32, 0.32),
  xlab = "Estimated interaction (truth is exactly 0)",
  main = "A  Guides that did nothing", cex.main = 1.1
)
abline(v = 0, lty = 2, lwd = 1.6, col = "#333333")

# ---- Panel B: what the bias costs ------------------------------------------
par(mar = c(4.6, 4.6, 3.0, 1))
x_at <- -log10(thresholds)
x_text <- as.expression(lapply(thresholds, function(value) {
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
axis(1, at = x_at, labels = x_text)
for (entry in series) {
  rows <- scan_summary[
    scan_summary$method == entry$key & scan_summary$metric == "f1", ,
    drop = FALSE
  ]
  rows <- rows[match(thresholds, rows$nominal_fdr), , drop = FALSE]
  lines(x_at, rows$mean, type = "o", pch = entry$pch, lwd = 2,
        col = entry$colour)
}

par(mar = c(0, 0, 0, 0))
plot.new()
legend(
  "center",
  legend = vapply(series, function(entry) entry$label, character(1)),
  col = vapply(series, function(entry) entry$colour, character(1)),
  pch = vapply(series, function(entry) entry$pch, numeric(1)),
  lty = 1, lwd = 2, ncol = 2, bty = "n", cex = 1.3
)
dev.off()

cat("Wrote figures/simcrispr_interaction.pdf\n")
