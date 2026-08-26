#!/usr/bin/env Rscript

# BARCS re-analysis of the Liu et al. in vivo T cell CRISPR screen
# (Nature 2026, doi:10.1038/s41586-026-10906-9, Supplementary Table 5).
#
# This is the runnable, BARCS-package-only replacement for the original
# BARCS Studio export. It reads the committed model inputs under
# `results/liu_tcell/input/` and regenerates every file under
# `results/liu_tcell/output/` with the BARCS R package directly, so the
# result tree can be rebuilt on a clean checkout without the Studio app.
#
# Pipeline, applied per arm and per pairing factor:
#   1. bb_screen(term = "gate")          guide-level beta-binomial regression
#   2. bb_calibrate_controls()           non-targeting-control tail calibration
#                                        (adds the raw_* columns to GUIDES)
#   3. bb_gene_consistency()             control-aware gene-level summary (GENES)
#
# The counts in Supplementary Table 5 are MAGeCK median-normalised and were
# already integer-valued and non-negative when exported, and each library's
# `total` in the metadata equals the column sum of its counts. They are used
# as the beta-binomial denominators unchanged. This is a re-analysis of
# published normalised counts, not an independent replication; see
# `results/liu_tcell/BARCS_summary.md` for the full interpretation and
# limits.

options(stringsAsFactors = FALSE)

source(file.path("R", "load_barcs.R"))

input_dir <- file.path("results", "liu_tcell", "input")
output_dir <- file.path("results", "liu_tcell", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Each arm is a separate MAGeCK normalisation run and is fitted on its own.
# `pairings` lists the design formulas evaluated for each arm: donor pairing is
# the primary model, and the Arm A mouse pairing is the Result 3 sensitivity
# fit. The `label` becomes the output file stem, matching the committed files.
arms <- list(
  armA = list(
    counts = "armA_counts.csv",
    metadata = "armA_metadata.csv",
    pairings = list(
      list(label = "armA_gate-donor", formula = ~ gate + donor),
      list(label = "armA_gate-mouse", formula = ~ gate + mouse)
    )
  ),
  armB = list(
    counts = "armB_counts.csv",
    metadata = "armB_metadata.csv",
    pairings = list(
      list(label = "armB_gate-donor", formula = ~ gate + donor)
    )
  )
)

read_arm <- function(arm) {
  counts_path <- file.path(input_dir, arm$counts)
  metadata_path <- file.path(input_dir, arm$metadata)
  if (!file.exists(counts_path) || !file.exists(metadata_path)) {
    stop(
      "Missing Liu T cell inputs: ", counts_path, " / ", metadata_path,
      call. = FALSE
    )
  }

  counts_table <- read.csv(counts_path, check.names = FALSE)
  metadata <- read.csv(metadata_path, check.names = FALSE)

  # `gate`, `donor`, and `mouse` are model factors, not continuous covariates.
  metadata$gate <- factor(metadata$gate)
  metadata$donor <- factor(metadata$donor)
  metadata$mouse <- factor(metadata$mouse)

  count_matrix <- as.matrix(counts_table[, metadata$sample, drop = FALSE])
  storage.mode(count_matrix) <- "double"

  # A beta-binomial denominator counts reads, so non-integer counts would make
  # every fit silently return NA; fail loudly if the inputs ever regress.
  if (any(abs(count_matrix - round(count_matrix)) >= sqrt(.Machine$double.eps))) {
    stop("Input counts must be integer-valued library counts.", call. = FALSE)
  }

  list(
    counts = count_matrix,
    metadata = metadata,
    guide = counts_table$guide,
    gene = counts_table$gene,
    control = tolower(as.character(counts_table$control)) == "true",
    totals = metadata$total
  )
}

fit_pairing <- function(arm_data, pairing) {
  guides <- bb_screen(
    counts = arm_data$counts,
    data = arm_data$metadata,
    formula = pairing$formula,
    term = "gate1",
    totals = arm_data$totals,
    guide = arm_data$guide,
    gene = arm_data$gene,
    min_total_count = 0
  )

  guides <- bb_calibrate_controls(guides, control = arm_data$control)

  guides$control <- arm_data$control
  guide_columns <- c(
    "gene", "guide", "estimate", "std_error", "t_value", "df", "p_value",
    "fdr", "rho", "pearson_null", "mean_cpm", "converged", "control",
    "raw_std_error", "raw_t_value", "raw_p_value", "raw_fdr"
  )
  guides <- guides[, guide_columns, drop = FALSE]

  genes <- bb_gene_consistency(guides, control = arm_data$control)

  # null_center / null_scale come back as attributes; promote them to columns
  # so the shared calibration constants are visible in the written table.
  genes$null_center <- attr(genes, "null_center")
  genes$null_scale <- attr(genes, "null_scale")

  gene_columns <- c(
    "gene", "n_guides", "estimate", "std_error", "statistic", "p_value",
    "fdr", "guide_direction_agreement", "raw_statistic", "converged_fraction",
    "control_gene", "null_center", "null_scale"
  )
  genes <- genes[, gene_columns, drop = FALSE]

  list(guides = guides, genes = genes)
}

for (arm_name in names(arms)) {
  arm <- arms[[arm_name]]
  arm_data <- read_arm(arm)

  for (pairing in arm$pairings) {
    fit <- fit_pairing(arm_data, pairing)

    guides_path <- file.path(output_dir, paste0(pairing$label, "_GUIDES.csv"))
    genes_path <- file.path(output_dir, paste0(pairing$label, "_GENES.csv"))
    write.csv(fit$guides, guides_path, row.names = FALSE)
    write.csv(fit$genes, genes_path, row.names = FALSE)

    message(sprintf(
      "%s: %d guides, %d genes -> %s, %s",
      pairing$label, nrow(fit$guides), nrow(fit$genes),
      basename(guides_path), basename(genes_path)
    ))
  }
}
