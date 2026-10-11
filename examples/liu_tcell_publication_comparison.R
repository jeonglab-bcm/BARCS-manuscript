#!/usr/bin/env Rscript

# BARCS versus the published Liu et al. analysis of the in vivo T cell IFNγ
# screen (Nature 2026, doi:10.1038/s41586-026-10906-9, Extended Data Fig. 5).
#
# The paper summarises the screen with MAGeCK: a one-sided RRA test and a
# per-gene log2 fold change between the IFNγ-high and IFNγ-low sort gates.
# BARCS instead reports a logit-scale beta-binomial coefficient with a
# calibrated p-value. The two are comparable in sign and rank, not in units,
# so every panel is read as concordance, never as a fitted slope.
#
# Panels, mapped to Extended Data Fig. 5:
#   A  BARCS effect vs MAGeCK log2FC, Arm A (CD3-scFv vs A375low) — method
#      concordance on the same screen the paper's panel f volcano summarises.
#   B  the same for Arm B (NY-ESO-1 TCR vs WT A375).
#   C  BARCS Arm A vs Arm B effect — the cross-model agreement the paper shows
#      with RRA scores in panel g, redrawn with BARCS.
#   D  BARCS Arm A volcano (effect vs -log10 p), the BARCS counterpart of the
#      MAGeCK RRA volcano in panel f.
#
# Requires results/liu_tcell/output/ (run examples/liu_tcell_barcs.R) and the
# published workbook under data/raw/liu_tcell/ (run scripts/prepare_liu_tcell.R).

options(stringsAsFactors = FALSE)

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("The `readxl` package is required.")
}

output_dir <- file.path("results", "liu_tcell", "output")
table5_path <- file.path("data", "raw", "liu_tcell", "Supplementary_Table_5.xlsx")
if (!file.exists(table5_path)) {
  stop(
    "Missing ", table5_path, ". Run `Rscript scripts/prepare_liu_tcell.R` first.",
    call. = FALSE
  )
}

read_barcs_genes <- function(stem) {
  path <- file.path(output_dir, paste0(stem, "_GENES.csv"))
  if (!file.exists(path)) {
    stop(
      "Missing ", path, ". Run `Rscript examples/liu_tcell_barcs.R` first.",
      call. = FALSE
    )
  }
  read.csv(path)
}

read_published_genes <- function(sheet) {
  published <- as.data.frame(
    suppressMessages(readxl::read_excel(table5_path, sheet = sheet)),
    check.names = FALSE
  )
  # MAGeCK reports the same per-gene log2FC in the neg| and pos| blocks; the
  # one-sided RRA p-value is taken from whichever tail the gene falls in.
  enriched <- published[["pos|p-value"]] <= published[["neg|p-value"]]
  data.frame(
    gene = published$id,
    log2fc = published[["pos|lfc"]],
    rra_p = ifelse(enriched, published[["pos|p-value"]], published[["neg|p-value"]]),
    rra_enriched = enriched
  )
}

barcs_armA_donor <- read_barcs_genes("armA_gate-donor")
barcs_armB_donor <- read_barcs_genes("armB_gate-donor")
published_armA <- read_published_genes("CD3scFv_A375low_3Donors_geneSum")
published_armB <- read_published_genes("WT_A375_2donors_geneSum")

# Okabe-Ito, matching R/method_palette.R: BARCS blue, MAGeCK orange, grey null.
colour_hit <- "#0072B2"
colour_enriched <- "#D55E00"
colour_depleted <- "#0072B2"
colour_null <- "#8C8C8C"
hit_fdr <- 0.10

# Genes worth naming: the discussed hits plus the two the pairing analysis
# turns on. Labelling every gene would bury the ones that matter.
label_genes <- c(
  "IFNG", "TNF", "GNAS", "STUB1", "TNFAIP3", "RASA2", "NFKBIA", "MED12",
  "STT3B", "P2RY8", "PTGER4"
)

# Spread labels around a dense point so names in a tight cluster do not
# overprint: pos cycles below / left / above / right by gene order.
label_positions <- function(genes) {
  positions <- c(1, 2, 3, 4)
  positions[(match(genes, sort(unique(genes))) - 1L) %% 4L + 1L]
}

# Pad a range so gene labels next to edge points have room and are not clipped.
padded_range <- function(values, fraction = 0.12) {
  span <- diff(range(values))
  range(values) + c(-1, 1) * fraction * span
}

concordance_panel <- function(barcs, published, letter, subtitle) {
  merged <- merge(
    barcs[, c("gene", "estimate", "fdr")], published[, c("gene", "log2fc")],
    by = "gene"
  )
  hit <- merged$fdr < hit_fdr
  spearman <- cor(merged$estimate, merged$log2fc, method = "spearman")
  agreement <- mean(sign(merged$estimate) == sign(merged$log2fc))

  par(mar = c(4.1, 4.2, 2.7, 1.4), cex.axis = 0.92, cex.lab = 0.98)
  plot(
    merged$estimate, merged$log2fc,
    type = "n",
    xlim = padded_range(merged$estimate),
    ylim = padded_range(merged$log2fc),
    xlab = "BARCS effect (logit-scale coefficient)",
    ylab = "MAGeCK gene log2 fold change",
    main = paste0(letter, "  ", subtitle),
    bty = "l", cex.main = 1.05
  )
  abline(h = 0, col = "#666666", lty = 3)
  abline(v = 0, col = "#666666", lty = 3)
  points(
    merged$estimate[!hit], merged$log2fc[!hit],
    pch = 16, cex = 0.9, col = adjustcolor(colour_null, alpha.f = 0.55)
  )
  points(
    merged$estimate[hit], merged$log2fc[hit],
    pch = 16, cex = 1.1, col = adjustcolor(colour_hit, alpha.f = 0.85)
  )
  to_label <- merged$gene %in% label_genes
  text(
    merged$estimate[to_label], merged$log2fc[to_label],
    labels = merged$gene[to_label],
    pos = label_positions(merged$gene[to_label]), offset = 0.4, cex = 0.7,
    col = "#333333", xpd = NA
  )
  legend(
    "topleft",
    legend = c(
      sprintf("Spearman rho = %.2f", spearman),
      sprintf("sign agreement = %.0f%%", 100 * agreement),
      "Called by BARCS (FDR < 0.10)"
    ),
    pch = c(NA, NA, 16), col = c(NA, NA, colour_hit),
    bty = "n", cex = 0.8, text.col = "#333333"
  )
}

cross_arm_panel <- function(armA, armB, letter) {
  merged <- merge(
    armA[, c("gene", "estimate", "fdr")],
    armB[, c("gene", "estimate", "fdr")],
    by = "gene", suffixes = c("_A", "_B")
  )
  merged <- merged[merged$gene != "NTCTRL", ]
  hit <- merged$fdr_A < hit_fdr | merged$fdr_B < hit_fdr
  pearson <- cor(merged$estimate_A, merged$estimate_B)

  par(mar = c(4.1, 4.2, 2.7, 1.4), cex.axis = 0.92, cex.lab = 0.98)
  plot(
    merged$estimate_A, merged$estimate_B,
    type = "n",
    xlim = padded_range(merged$estimate_A),
    ylim = padded_range(merged$estimate_B),
    xlab = "BARCS effect - Arm A (A375low)",
    ylab = "BARCS effect - Arm B (NY-ESO-1)",
    main = paste0(letter, "  Cross-model agreement (cf. Fig. 5g)"),
    bty = "l", cex.main = 1.05
  )
  abline(h = 0, col = "#666666", lty = 3)
  abline(v = 0, col = "#666666", lty = 3)
  points(
    merged$estimate_A[!hit], merged$estimate_B[!hit],
    pch = 16, cex = 0.9, col = adjustcolor(colour_null, alpha.f = 0.55)
  )
  points(
    merged$estimate_A[hit], merged$estimate_B[hit],
    pch = 16, cex = 1.1, col = adjustcolor(colour_hit, alpha.f = 0.85)
  )
  to_label <- merged$gene %in% label_genes
  text(
    merged$estimate_A[to_label], merged$estimate_B[to_label],
    labels = merged$gene[to_label],
    pos = label_positions(merged$gene[to_label]), offset = 0.4, cex = 0.7,
    col = "#333333", xpd = NA
  )
  legend(
    "topleft",
    legend = c(
      sprintf("Pearson r = %.2f", pearson),
      "Called by BARCS (FDR < 0.10)"
    ),
    pch = c(NA, 16), col = c(NA, colour_hit),
    bty = "n", cex = 0.8, text.col = "#333333"
  )
}

volcano_panel <- function(barcs, letter) {
  effect <- barcs$estimate
  y <- -log10(barcs$p_value)
  hit <- barcs$fdr < hit_fdr
  base_colour <- ifelse(
    !hit, colour_null,
    ifelse(effect > 0, colour_enriched, colour_depleted)
  )
  point_colour <- mapply(
    function(col, alpha) adjustcolor(col, alpha.f = alpha),
    base_colour, ifelse(hit, 0.85, 0.5)
  )

  par(mar = c(4.1, 4.2, 2.7, 1.4), cex.axis = 0.92, cex.lab = 0.98)
  plot(
    effect, y,
    type = "n",
    xlim = padded_range(effect),
    ylim = c(0, max(y) * 1.05),
    xlab = "BARCS effect (logit-scale coefficient)",
    ylab = expression(-log[10](italic(p))),
    main = paste0(letter, "  BARCS Arm A volcano (cf. RRA Fig. 5f)"),
    bty = "l", cex.main = 1.05
  )
  abline(v = 0, col = "#666666", lty = 3)
  abline(h = -log10(0.05), col = "#3B5BDB", lty = 2)
  points(
    effect, y, pch = 16,
    cex = ifelse(hit, 1.1, 0.9),
    col = point_colour
  )
  to_label <- barcs$gene %in% label_genes & hit
  text(
    effect[to_label], y[to_label],
    labels = barcs$gene[to_label],
    pos = label_positions(barcs$gene[to_label]), offset = 0.4, cex = 0.7,
    col = "#333333", xpd = NA
  )
  legend(
    "topright",
    legend = c("Enriched (FDR < 0.10)", "Depleted (FDR < 0.10)", "Not called"),
    pch = 16, pt.cex = 1.0, cex = 0.72, bty = "n", text.col = "#333333",
    col = c(colour_enriched, colour_depleted, colour_null)
  )
}

draw_figure <- function() {
  layout(matrix(seq_len(4), nrow = 2, byrow = TRUE))
  concordance_panel(
    barcs_armA_donor, published_armA, "a", "Arm A: BARCS vs MAGeCK (A375low)"
  )
  concordance_panel(
    barcs_armB_donor, published_armB, "b", "Arm B: BARCS vs MAGeCK (NY-ESO-1)"
  )
  cross_arm_panel(barcs_armA_donor, barcs_armB_donor, "c")
  volcano_panel(barcs_armA_donor, "d")
}

dir.create("figures", showWarnings = FALSE)
figure_stem <- file.path("figures", "liu_tcell_publication_comparison")

# PNG for quick review (a repo-only diagnostic), PDF for the manuscript, which
# includes figures as vector PDF.
png(paste0(figure_stem, ".png"), width = 2000, height = 2000, res = 200)
draw_figure()
invisible(dev.off())

pdf(paste0(figure_stem, ".pdf"), width = 10, height = 10, pointsize = 12,
    useDingbats = FALSE)
draw_figure()
invisible(dev.off())

message("Wrote ", figure_stem, ".png and ", figure_stem, ".pdf")
