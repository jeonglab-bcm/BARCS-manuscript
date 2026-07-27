#!/usr/bin/env Rscript

# Reference fits for the Python port's equivalence tests.
#
# `python/tests/test_against_r.py` writes a fixture screen to a directory,
# runs this script over it, and compares the CSVs it writes against the same
# quantities computed by `python/barcs/`. The Python side owns the fixture so
# both implementations are handed byte-identical input; this script only
# fits and dumps.
#
#     Rscript python/tests/reference_fit.R <fixture_dir>
#
# Anything printed here is diagnostic only. The comparison reads the CSVs.

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) {
  stop("usage: reference_fit.R <fixture_dir>", call. = FALSE)
}
fixture_dir <- arguments[[1L]]

# Run from the repository root so this relative source path resolves.
source(file.path("R", "bbreg.R"))

counts <- as.matrix(read.csv(
  file.path(fixture_dir, "counts.csv"), row.names = 1L, check.names = FALSE
))
storage.mode(counts) <- "double"
samples <- read.csv(file.path(fixture_dir, "samples.csv"))
samples$batch <- factor(samples$batch)
totals <- scan(file.path(fixture_dir, "totals.txt"), quiet = TRUE)
guide_gene <- read.csv(file.path(fixture_dir, "guide_gene.csv"))

# 1. A single-guide fit, reported in full. This is the tightest comparison:
#    every coefficient, not just the reported term.
single <- bbreg(
  count = counts[1L, ],
  total = totals,
  formula = ~ dose + batch,
  data = samples
)
single_table <- as.data.frame(single$coefficient_table)
single_table$coefficient <- rownames(single_table)
single_table$rho <- single$rho
single_table$pearson <- single$pearson
single_table$pearson_null <- single$pearson_null
single_table$scale <- single$scale
single_table$converged <- single$converged
write.csv(
  single_table[, c("coefficient", "estimate", "std_error", "t_value", "df",
                   "p_value", "rho", "pearson", "pearson_null", "scale",
                   "converged")],
  file.path(fixture_dir, "r_single_fit.csv"),
  row.names = FALSE
)

# 2. A named contrast on that same fit: two dose steps.
contrast <- bb_contrast(single, c(dose = 2))
write.csv(contrast, file.path(fixture_dir, "r_contrast.csv"), row.names = FALSE)

# 3. The whole screen, guide by guide.
screen <- bb_screen(
  counts = counts,
  data = samples,
  formula = ~ dose + batch,
  term = "dose",
  totals = totals,
  guide = guide_gene$guide,
  gene = guide_gene$gene
)
write.csv(screen, file.path(fixture_dir, "r_screen.csv"), row.names = FALSE)

# 4. Negative-control calibration, both estimators.
control <- guide_gene$gene == "control"
for (method in c("tail_quantile", "qq_slope")) {
  calibrated <- bb_calibrate_controls(screen, control = control, method = method)
  out <- calibrated[, c("guide", "std_error", "t_value", "p_value", "fdr")]
  out$control_scale <- attr(calibrated, "control_scale")
  write.csv(
    out,
    file.path(fixture_dir, paste0("r_calibrated_", method, ".csv")),
    row.names = FALSE
  )
}

# 5. Guide-to-gene aggregation.
gene_result <- bb_gene_original(screen)
write.csv(
  gene_result, file.path(fixture_dir, "r_gene.csv"), row.names = FALSE
)

cat("reference fits written to", fixture_dir, "\n")
