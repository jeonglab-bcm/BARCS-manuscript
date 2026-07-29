#!/usr/bin/env Rscript

# Minimal BARCS analysis.
#
# This example is deliberately small enough to read in one pass. It shows the
# three objects that every analysis needs:
#   1. a guide-by-library count matrix;
#   2. immutable full-library totals; and
#   3. one row of sample metadata per count-matrix column.

if (requireNamespace("BARCS", quietly = TRUE)) {
  suppressPackageStartupMessages(library(BARCS))
} else {
  source(file.path("R", "bbreg.R"))
}

counts <- rbind(
  guide_1 = c(120, 110, 85, 70, 48, 33),
  guide_2 = c(75,  79,  72, 70, 68, 66),
  control = c(42,  45,  41, 43, 44, 40)
)
colnames(counts) <- paste0("sample_", seq_len(ncol(counts)))

# These totals represent all mapped guides in each library, not only the three
# rows retained in this toy count matrix.
library_totals <- c(100000, 98000, 105000, 101000, 99000, 103000)

sample_data <- data.frame(
  sample = colnames(counts),
  replicate = factor(rep(c("A", "B"), times = 3)),
  day = rep(c(0, 7, 14), each = 2)
)

result <- bb_screen(
  counts = counts,
  totals = library_totals,
  data = sample_data,
  formula = ~ replicate + I(day / 14),
  term = "I(day/14)",
  guide = rownames(counts),
  min_total_count = 0
)

print(result[, c("guide", "estimate", "std_error", "p_value", "fdr")])
