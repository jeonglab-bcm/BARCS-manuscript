#!/usr/bin/env Rscript

# Held-out null calibration for BARCS and the published MAGeCK pipeline on the
# Schmidt et al. T cell screens (Science 2022, doi:10.1126/science.abj4008).
#
# On this screen the two methods recover the published biology about equally
# well but disagree by a factor of three on how many genes they call. Recovery
# cannot say which list is right, because it never asks what either method does
# with a gene that has no effect. This does.
#
# A screen's 992 nontargeting sgRNAs are grouped into 165 gene-shaped
# pseudo-genes (R/schmidt_tcell.R). Alternating pseudo-genes form a calibration
# half and an evaluation half, and the two folds swap those roles. Within a
# fold, the calibration half is the only control set either method may tune its
# null on -- BARCS's control calibration, MAGeCK's --control-sgrna -- and the
# evaluation half is relabelled and scored alongside the ~18,800 real genes.
# Every evaluation pseudo-gene is null by construction, so anything called among
# them is a false positive; across both folds each is evaluated exactly once
# while held out. Both methods see the identical grouping and the identical
# folds, and neither is scored on controls it was calibrated on.
#
# A third arm justifies the dispersion-moderation step in
# examples/schmidt_tcell_barcs.R. `bb_moderate_dispersion()` is not free to
# assume -- on the 268-guide Liu screen the same step costs published hits --
# so the "BARCS unmoderated" arm rescores the identical fits from the
# pre-moderation columns the guide table preserves, and the held-out controls
# say whether moderation bought its extra discoveries out of the null.
#
# Requires results/schmidt_tcell/output/ (run examples/schmidt_tcell_barcs.R)
# and the `mageck` executable; set MAGECK to its path if it is not on PATH.
# Writes results/schmidt_tcell/comparison/null_calibration.csv.

options(stringsAsFactors = FALSE)

source(file.path("R", "load_barcs.R"))
source(file.path("R", "schmidt_tcell.R"))

input_dir <- file.path("results", "schmidt_tcell", "input")
output_dir <- file.path("results", "schmidt_tcell", "output")
comparison_dir <- file.path("results", "schmidt_tcell", "comparison")
holdout_dir <- file.path("results", "schmidt_tcell", "mageck_holdout")
dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)

folds <- c(0L, 1L)

# Relabel a screen's control guides for one fold. The calibration half keeps
# control status; the evaluation half loses it and becomes ordinary-looking
# genes, so every method scores it exactly as it scores a real gene.
relabel_controls <- function(gene, guide, is_control, fold) {
  index <- schmidt_control_pseudogene(guide[is_control])
  role <- schmidt_control_fold(index, fold)
  gene[is_control] <- sprintf(
    "%s_%03d", ifelse(role == "calibration", "NTCTRL", "NTEVAL"), index
  )
  list(
    gene = gene,
    control = replace(is_control, is_control, role == "calibration")
  )
}

rows <- list()
for (screen in schmidt_screens) {
  parts <- strsplit(screen, "_", fixed = TRUE)[[1L]]

  # The guide-level beta-binomial fits do not depend on gene labels, so the
  # committed guide table is reused; only the control calibration and the
  # gene-level summary are redone per fold. The calibrated columns are rolled
  # back to the raw ones first, so each fold's calibration is applied once, to
  # that fold's calibration half alone.
  fitted <- read.csv(file.path(output_dir, paste0(screen, "_GUIDES.csv.gz")))
  # Roll the control calibration back so each fold applies its own, once.
  fitted$std_error <- fitted$raw_std_error
  fitted$t_value <- fitted$raw_t_value
  fitted$p_value <- fitted$raw_p_value
  fitted$fdr <- fitted$raw_fdr
  fitted[c("raw_std_error", "raw_t_value", "raw_p_value", "raw_fdr")] <- NULL

  # The pre-moderation statistics `bb_moderate_dispersion()` preserved, which
  # give the unmoderated arm without refitting.
  unmoderated <- fitted
  unmoderated$std_error <- fitted$unmoderated_std_error
  unmoderated$t_value <- fitted$unmoderated_t_value
  unmoderated$df <- fitted$unmoderated_df
  unmoderated$p_value <- fitted$unmoderated_p_value
  unmoderated$fdr <- p.adjust(unmoderated$p_value, method = "BH")

  published <- schmidt_mageck_table(input_dir, parts[1L], parts[2L])

  for (fold in folds) {
    # ---- BARCS, moderated (the shipped pipeline) and unmoderated -----------
    score_barcs <- function(guides) {
      relabelled <- relabel_controls(
        guides$gene, guides$guide, guides$control, fold
      )
      guides$gene <- relabelled$gene
      guides$control <- relabelled$control
      guides <- bb_calibrate_controls(guides, control = guides$control)
      bb_gene_consistency(guides, control = guides$control)
    }
    barcs_genes <- score_barcs(fitted)
    barcs_eval <- barcs_genes[startsWith(barcs_genes$gene, "NTEVAL_"), ]
    plain_genes <- score_barcs(unmoderated)
    plain_eval <- plain_genes[startsWith(plain_genes$gene, "NTEVAL_"), ]

    # ---- MAGeCK ------------------------------------------------------------
    table <- published
    relabelled <- relabel_controls(
      table$Gene, table$sgRNA, table$Gene == "NO-TARGET", fold
    )
    table$Gene <- relabelled$gene

    mageck_genes <- schmidt_run_mageck(
      table = table,
      controls = table$sgRNA[relabelled$control],
      screen_dir = file.path(holdout_dir, screen),
      prefix = sprintf("%s_fold%d", screen, fold)
    )
    mageck_genes$p_value <- pmin(
      mageck_genes[["pos|p-value"]], mageck_genes[["neg|p-value"]]
    )
    mageck_genes$hit <- schmidt_mageck_hits(mageck_genes)
    mageck_eval <- mageck_genes[startsWith(mageck_genes$id, "NTEVAL_"), ]

    rows[[length(rows) + 1L]] <- rbind(
      data.frame(
        screen = screen, fold = fold, method = "BARCS",
        null_scale = attr(barcs_genes, "null_scale"),
        global_scale = attr(barcs_genes, "global_scale"),
        control_scale = attr(barcs_genes, "control_scale"),
        genes_called = sum(barcs_genes$fdr < 0.05, na.rm = TRUE),
        null_genes = nrow(barcs_eval),
        null_called = sum(barcs_eval$fdr < 0.05, na.rm = TRUE),
        null_p05 = sum(barcs_eval$p_value < 0.05, na.rm = TRUE),
        # BARCS reports one two-sided p-value, so 5% of null genes should fall
        # below 0.05. MAGeCK reports two one-sided RRA p-values and the
        # published rule reads whichever is smaller, and for a null gene the two
        # tails are near-complementary, so its matched expectation is 10%.
        null_p05_expected = 0.05
      ),
      data.frame(
        screen = screen, fold = fold, method = "BARCS unmoderated",
        null_scale = attr(plain_genes, "null_scale"),
        global_scale = attr(plain_genes, "global_scale"),
        control_scale = attr(plain_genes, "control_scale"),
        genes_called = sum(plain_genes$fdr < 0.05, na.rm = TRUE),
        null_genes = nrow(plain_eval),
        null_called = sum(plain_eval$fdr < 0.05, na.rm = TRUE),
        null_p05 = sum(plain_eval$p_value < 0.05, na.rm = TRUE),
        null_p05_expected = 0.05
      ),
      data.frame(
        screen = screen, fold = fold, method = "MAGeCK",
        null_scale = NA_real_, global_scale = NA_real_, control_scale = NA_real_,
        genes_called = sum(mageck_genes$hit, na.rm = TRUE),
        null_genes = nrow(mageck_eval),
        null_called = sum(mageck_eval$hit, na.rm = TRUE),
        null_p05 = sum(mageck_eval$p_value < 0.05, na.rm = TRUE),
        null_p05_expected = 0.10
      )
    )
  }
}

calibration <- do.call(rbind, rows)
rownames(calibration) <- NULL
calibration$null_p05_rate <- calibration$null_p05 / calibration$null_genes
calibration$null_p05_ratio <- calibration$null_p05_rate / calibration$null_p05_expected
write.csv(calibration, file.path(comparison_dir, "null_calibration.csv"),
          row.names = FALSE)

print(calibration[c("screen", "fold", "method", "null_scale", "genes_called",
                    "null_genes", "null_called", "null_p05", "null_p05_rate",
                    "null_p05_ratio")], row.names = FALSE)

# How much the hit count moves when only the held-out half of the controls is
# available to calibrate on. The control-derived scale is a tail quantile of at
# most 165 control-gene statistics, so it is a noisy estimate, and a null scale
# multiplies every gene statistic in the screen.
cat("\nHit count across the eight held-out fits, by method:\n")
print(aggregate(genes_called ~ method, calibration,
                function(v) c(min = min(v), median = median(v), max = max(v))))

pooled <- aggregate(
  cbind(genes_called, null_genes, null_called, null_p05) ~ method,
  data = calibration, FUN = sum
)
pooled$null_p05_rate <- pooled$null_p05 / pooled$null_genes
pooled$null_p05_ratio <- pooled$null_p05_rate /
  c(BARCS = 0.05, `BARCS unmoderated` = 0.05, MAGeCK = 0.10)[pooled$method]
cat("\nPooled over four screens and both folds:\n")
print(pooled, row.names = FALSE)
