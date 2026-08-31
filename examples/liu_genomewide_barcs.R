#!/usr/bin/env Rscript

# BARCS re-analysis of the genome-wide Liu et al. in vivo T cell IFNγ screen
# (Nature 2026, doi:10.1038/s41586-026-10906-9, Supplementary Table 2, the
# screen behind main Fig. 4).
#
# Reads the committed model inputs under results/liu_genomewide/input/ and
# writes gene- and guide-level BARCS results under results/liu_genomewide/output/.
# The design is `~ gate + donor`: the sort gate (IFNγ-high vs IFNγ-low) is the
# tested term and donor is kept as a covariate. Unlike the focused sub-library
# (examples/liu_tcell_barcs.R), the genome-wide screen is collapsed to two donor
# replicates, so there is no per-mouse pairing to fit here.
#
# Pipeline: bb_screen(term = "gate") -> bb_calibrate_controls() against the
# non-targeting (NTCTRL) guides -> bb_gene_consistency().
#
# Requires results/liu_genomewide/input/ (run scripts/prepare_liu_genomewide.R).
# This is a re-analysis of published median-normalised counts, not a
# replication; see results/liu_genomewide/README.md.

options(stringsAsFactors = FALSE)

source(file.path("R", "load_barcs.R"))

input_dir <- file.path("results", "liu_genomewide", "input")
output_dir <- file.path("results", "liu_genomewide", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Fitting ~77k guides is minutes on one core; use the spare cores if present.
ncores <- max(1L, min(8L, parallel::detectCores() - 2L))

counts_path <- file.path(input_dir, "counts.csv")
metadata_path <- file.path(input_dir, "metadata.csv")
if (!file.exists(counts_path) || !file.exists(metadata_path)) {
  stop(
    "Missing genome-wide inputs. Run `Rscript scripts/prepare_liu_genomewide.R` first.",
    call. = FALSE
  )
}

counts_table <- read.csv(counts_path, check.names = FALSE)
metadata <- read.csv(metadata_path, check.names = FALSE)
metadata$gate <- factor(metadata$gate)
metadata$donor <- factor(metadata$donor)

count_matrix <- as.matrix(counts_table[, metadata$sample, drop = FALSE])
storage.mode(count_matrix) <- "double"
if (any(abs(count_matrix - round(count_matrix)) >= sqrt(.Machine$double.eps))) {
  stop("Input counts must be integer-valued library counts.", call. = FALSE)
}

control <- tolower(as.character(counts_table$control)) == "true"

guides <- bb_screen(
  counts = count_matrix,
  data = metadata,
  formula = ~ gate + donor,
  term = "gate1",
  totals = metadata$total,
  guide = counts_table$guide,
  gene = counts_table$gene,
  min_total_count = 0,
  ncores = ncores
)

guides <- bb_calibrate_controls(guides, control = control)

guides$control <- control
guide_columns <- c(
  "gene", "guide", "estimate", "std_error", "t_value", "df", "p_value",
  "fdr", "rho", "pearson_null", "mean_cpm", "converged", "control",
  "raw_std_error", "raw_t_value", "raw_p_value", "raw_fdr"
)
guides <- guides[, guide_columns, drop = FALSE]

genes <- bb_gene_consistency(guides, control = control)

# null_center / null_scale come back as attributes; promote them to columns so
# the shared calibration constants are visible in the written table.
genes$null_center <- attr(genes, "null_center")
genes$null_scale <- attr(genes, "null_scale")
gene_columns <- c(
  "gene", "n_guides", "estimate", "std_error", "statistic", "p_value",
  "fdr", "guide_direction_agreement", "raw_statistic", "converged_fraction",
  "control_gene", "null_center", "null_scale"
)
genes <- genes[, gene_columns, drop = FALSE]

guides_path <- file.path(output_dir, "gate-donor_GUIDES.csv")
genes_path <- file.path(output_dir, "gate-donor_GENES.csv")
write.csv(guides, guides_path, row.names = FALSE)
write.csv(genes, genes_path, row.names = FALSE)

message(sprintf(
  "genome-wide gate-donor: %d guides, %d genes -> %s, %s",
  nrow(guides), nrow(genes), basename(guides_path), basename(genes_path)
))
message(sprintf(
  "convergence %.1f%%, control-calibration scale %.3f",
  100 * mean(guides$converged, na.rm = TRUE),
  attr(guides, "control_scale")
))
