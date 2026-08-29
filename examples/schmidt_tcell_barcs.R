#!/usr/bin/env Rscript

# BARCS re-analysis of the Schmidt et al. genome-wide CRISPRa/CRISPRi screens
# in primary human T cells (Science 2022, doi:10.1126/science.abj4008; raw
# counts GEO GSE174255).
#
# Reads the per-library-set inputs written by scripts/prepare_schmidt_tcell.R
# and writes gene- and guide-level results for the four screens
# (CRISPRa/CRISPRi x IL2/IFNG) into results/schmidt_tcell/output/.
#
# The design point of the re-analysis
# -----------------------------------
# The published analysis normalised every library to its total read count,
# merged the two library sets into one table, and ran `mageck test` once per
# screen with --paired. That pipeline has to put two separately sequenced pools
# on a common scale before it can compare anything.
#
# BARCS never needs to. A guide's high-versus-low effect is a contrast internal
# to the libraries that guide was sequenced in, against those libraries' own
# read totals, so guides from Set A and Set B are directly comparable without
# either set being rescaled onto the other. Each library set is therefore fitted
# on its own raw counts, and the two sets are only combined at the gene level,
# where they supply 3 + 3 = 6 guides per gene.
#
# Within a set, one model per guide covers all twelve sorted and unsorted
# libraries:
#
#   ~ donor + assay * bin      bin in {low, unsorted, high}, assay in {IFNG, IL2}
#
# so donor is an explicit term, the unsorted bins anchor the dispersion, and
# both cytokine contrasts are read off a single fit with 5 residual degrees of
# freedom. `bb_screen()` reports one model-matrix coefficient, so the fit is run
# once per cytokine with `assay` releveled; the two runs are the same model in
# two parameterisations, and `binhigh` is that cytokine's high-versus-low
# contrast each time.
#
# A positive effect means the perturbation enriched the guide in the
# cytokine-high bin, matching the published log2(high/low) sign convention.

options(stringsAsFactors = FALSE)

source(file.path("R", "load_barcs.R"))
source(file.path("R", "schmidt_tcell.R"))

input_dir <- file.path("results", "schmidt_tcell", "input")
output_dir <- file.path("results", "schmidt_tcell", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ncores <- max(1L, as.integer(Sys.getenv("BARCS_NCORES", unset = "4")))

assays <- c("IL2", "IFNG")

# One fit per (library set, cytokine): `assay` is releveled so that the cytokine
# being tested is the reference level, which makes `binhigh` its own
# high-versus-low contrast.
fit_set_assay <- function(set_data, assay) {
  data <- set_data$metadata
  data$assay <- relevel(factor(data$assay), ref = assay)

  guides <- bb_screen(
    counts = set_data$counts,
    data = data,
    formula = ~ donor + assay * bin,
    term = "binhigh",
    totals = set_data$totals,
    guide = set_data$guide,
    gene = set_data$gene,
    ncores = ncores
  )
  guides$control <- set_data$control
  guides$library_set <- set_data$set_name
  # Qualify the guide id with its library set: the two sets share no sequence,
  # but which pool a guide was sequenced in is what its denominators come from
  # and is what the pseudo-gene grouping below is computed within.
  guides$guide <- paste0(guides$guide, "|", sub("^.*_", "", set_data$set_name))
  # Controls are grouped into gene-shaped pseudo-genes (see R/schmidt_tcell.R):
  # a single 992-guide NO-TARGET gene cannot report how often a gene-sized group
  # of null guides gets called, and it leaves `bb_gene_consistency()` with one
  # control gene, too few to estimate a control-based null at all.
  is_control <- guides$control
  guides$gene[is_control] <- sprintf(
    "NO-TARGET_%03d", schmidt_control_pseudogene(guides$guide[is_control])
  )
  guides
}

guide_columns <- c(
  "library_set", "gene", "guide", "estimate", "std_error", "t_value", "df",
  "p_value", "fdr", "rho", "pearson_null", "mean_cpm", "converged", "control",
  "raw_std_error", "raw_t_value", "raw_p_value", "raw_fdr"
)
# Columns `bb_moderate_dispersion()` adds. It keeps the pre-moderation
# statistics, so carrying these through makes the step auditable and lets
# examples/schmidt_tcell_null_calibration.R score the unmoderated pipeline from
# the same written table without refitting.
moderation_columns <- c(
  "unmoderated_std_error", "unmoderated_t_value", "unmoderated_df",
  "unmoderated_p_value", "variance_inflation", "moderated_inflation"
)
gene_columns <- c(
  "gene", "n_guides", "estimate", "std_error", "statistic", "p_value", "fdr",
  "guide_direction_agreement", "raw_statistic", "converged_fraction",
  "control_gene", "null_center", "null_scale", "global_scale", "control_scale"
)

set_fits <- list()
for (set_name in schmidt_library_sets) {
  set_data <- schmidt_read_set(input_dir, set_name)
  for (assay in assays) {
    started <- Sys.time()
    set_fits[[paste(set_name, assay, sep = "|")]] <- fit_set_assay(set_data, assay)
    message(sprintf(
      "fitted %s %s: %d guides in %.0f s",
      set_name, assay, nrow(set_data$counts),
      as.numeric(difftime(Sys.time(), started, units = "secs"))
    ))
  }
  rm(set_data)
}

for (modality in c("CRISPRa", "CRISPRi")) {
  for (assay in assays) {
    # Stack the two library sets: 3 guides per gene from each, on their own
    # denominators, never rescaled onto a shared normalisation.
    guides <- do.call(rbind, lapply(
      paste(paste0(modality, c("_SetA", "_SetB")), assay, sep = "|"),
      function(key) set_fits[[key]]
    ))

    # Empirical-Bayes moderation of the per-guide dispersions, borrowing across
    # the ~113,000 guides that share this design. At 5 residual degrees of
    # freedom a per-guide dispersion is badly determined, and leaving it
    # unmoderated is a real cost: on this screen moderation raises held-out
    # discoveries from 1,677 to 2,082 while *lowering* the held-out null rate
    # from 0.036 to 0.026, so it buys power and calibration together rather
    # than trading one for the other (see the null-calibration script).
    #
    # This is deliberately not treated as a universal default. On the Liu
    # screen in results/liu_tcell/ -- 268 guides, 36 genes -- the same step
    # costs four published hits and raises the null scale, because the prior is
    # estimated from too few guides. It is applied here because this screen has
    # the guide count to support it.
    guides <- bb_moderate_dispersion(guides)
    guides <- bb_calibrate_controls(guides, control = guides$control)
    guides <- guides[, c(guide_columns,
                         intersect(moderation_columns, names(guides))),
                     drop = FALSE]

    genes <- bb_gene_consistency(guides, control = guides$control)
    genes$null_center <- attr(genes, "null_center")
    genes$null_scale <- attr(genes, "null_scale")
    # The two candidates `bb_gene_consistency()` takes the maximum of. Whichever
    # is larger sets the scale for the whole screen, and the control-derived one
    # is a tail quantile of 165 pseudo-genes, so it is the noisier of the two;
    # carrying both makes that visible in the written table.
    genes$global_scale <- attr(genes, "global_scale")
    genes$control_scale <- attr(genes, "control_scale")
    genes <- genes[, gene_columns, drop = FALSE]
    genes <- genes[order(genes$p_value, -abs(genes$estimate)), ]

    label <- paste(modality, assay, sep = "_")
    genes_path <- file.path(output_dir, paste0(label, "_GENES.csv.gz"))
    guides_path <- file.path(output_dir, paste0(label, "_GUIDES.csv.gz"))
    write.csv(genes, gzfile(genes_path), row.names = FALSE)
    write.csv(guides, gzfile(guides_path), row.names = FALSE)

    control_genes <- genes[genes$control_gene, ]
    message(sprintf(
      paste("%s: %d guides, %d genes, %d at FDR<0.05;",
            "null scale %.2f, %d nontargeting pseudo-genes, %d of them called"),
      label, nrow(guides), nrow(genes), sum(genes$fdr < 0.05, na.rm = TRUE),
      genes$null_scale[1L], nrow(control_genes),
      sum(control_genes$fdr < 0.05, na.rm = TRUE)
    ))
  }
}
