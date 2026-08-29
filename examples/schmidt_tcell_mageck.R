#!/usr/bin/env Rscript

# Reproduce the published MAGeCK analysis of the Schmidt et al. T cell screens
# (Science 2022, doi:10.1126/science.abj4008) on the same GEO raw counts BARCS
# is fitted to, so the two methods are compared on identical inputs.
#
# The Methods recipe, followed step for step:
#
#   "raw read counts across both library sets were normalized to the total read
#    count in each sample, and each of the matching samples across two sets were
#    merged to generate a single normalized read count table. Normalized read
#    counts in high versus low bins were compared using mageck test with
#    --norm-method none, --paired, and --control-sgrna options, pairing samples
#    by donor and using nontargeting sgRNAs as controls"
#
# The table construction and the `mageck test` call both live in
# R/schmidt_tcell.R, so the null-calibration check runs MAGeCK exactly the same
# way. The one detail the Methods leave implicit -- what "normalized to
# the total read count" scales to -- is documented there.
#
# Requires the `mageck` executable; set MAGECK to its path if it is not on PATH.
# Writes results/schmidt_tcell/mageck/<screen>/<screen>.gene_summary.txt.

options(stringsAsFactors = FALSE)

source(file.path("R", "schmidt_tcell.R"))

input_dir <- file.path("results", "schmidt_tcell", "input")
mageck_dir <- file.path("results", "schmidt_tcell", "mageck")
dir.create(mageck_dir, recursive = TRUE, showWarnings = FALSE)

for (screen in schmidt_screens) {
  parts <- strsplit(screen, "_", fixed = TRUE)[[1L]]
  table <- schmidt_mageck_table(input_dir, parts[1L], parts[2L])

  genes <- schmidt_run_mageck(
    table = table,
    controls = table$sgRNA[table$Gene == "NO-TARGET"],
    screen_dir = file.path(mageck_dir, screen),
    prefix = screen
  )

  message(sprintf(
    "%s: %d guides merged, %d genes, %d hits (|median lfc| > 0.5 & FDR < 0.05)",
    screen, nrow(table), nrow(genes), sum(schmidt_mageck_hits(genes), na.rm = TRUE)
  ))
}
