#!/usr/bin/env Rscript

# Exploratory diagnostic: BARCS versus published MAGeCK RRA significance,
# head to head, for the Liu et al. in vivo T cell IFNγ screen (Nature 2026,
# doi:10.1038/s41586-026-10906-9). Arm A (CD3-scFv vs A375low).
#
# This plots -log10 p and -log10 FDR for the two methods against each other.
# It is deliberately NOT part of the manuscript figure
# (examples/liu_tcell_publication_comparison.R) and should not be read as a
# method-agreement metric. MAGeCK's one-sided RRA p-values are computed by
# permutation and floor at ~5e-6: on Arm A, 23 of 36 genes sit exactly on that
# floor, so the significance axis is heavily tied. The tie block smears across
# the whole BARCS range at one fixed y, collapsing the rank correlation
# (Spearman ~0.4 on p, ~0.3 on FDR) even though the two callers agree well on
# effect size (Spearman ~0.90 on the BARCS effect vs MAGeCK log2FC comparison
# in the manuscript figure). The honest concordance comparison is that
# effect/log2FC one; this diagnostic exists only to make the RRA discretisation
# visible.
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

barcs_path <- file.path(output_dir, "armA_gate-donor_GENES.csv")
if (!file.exists(barcs_path)) {
  stop(
    "Missing ", barcs_path, ". Run `Rscript examples/liu_tcell_barcs.R` first.",
    call. = FALSE
  )
}

published <- as.data.frame(
  suppressMessages(readxl::read_excel(table5_path, sheet = "CD3scFv_A375low_3Donors_geneSum")),
  check.names = FALSE
)
# MAGeCK is one-sided in each direction; take the significant tail per gene.
published$rra_p <- pmin(published[["pos|p-value"]], published[["neg|p-value"]])
published$rra_fdr <- pmin(published[["pos|fdr"]], published[["neg|fdr"]])

barcs <- read.csv(barcs_path)

merged <- merge(
  barcs[, c("gene", "p_value", "fdr")],
  data.frame(
    gene = published$id, rra_p = published$rra_p, rra_fdr = published$rra_fdr
  ),
  by = "gene"
)
merged <- merged[merged$gene != "NTCTRL", ]

# Okabe-Ito, matching R/method_palette.R.
colour_hit <- "#0072B2"
colour_null <- "#8C8C8C"
colour_floor <- "#C0392B"
hit_fdr <- 0.10

label_genes <- c(
  "IFNG", "TNF", "GNAS", "STUB1", "TNFAIP3", "RASA2", "NFKBIA", "MED12",
  "STT3B", "P2RY8", "PTGER4"
)

headtohead_panel <- function(barcs_neglog, mageck_neglog, xlab, ylab, letter,
                             subtitle, floor_value) {
  hit <- merged$fdr < hit_fdr
  spearman <- cor(barcs_neglog, mageck_neglog, method = "spearman")
  at_floor <- sum(abs(mageck_neglog - floor_value) < 1e-9)

  par(mar = c(4.2, 4.4, 2.8, 1.2), cex.axis = 0.92, cex.lab = 0.98)
  plot(
    barcs_neglog, mageck_neglog, type = "n",
    xlim = c(0, max(barcs_neglog) * 1.08),
    ylim = c(0, max(mageck_neglog) * 1.12),
    xlab = xlab, ylab = ylab,
    main = paste0(letter, "  ", subtitle),
    bty = "l", cex.main = 1.02
  )
  # The MAGeCK RRA permutation resolution limit, where the ties accumulate.
  abline(h = floor_value, col = colour_floor, lty = 2)
  text(
    max(barcs_neglog) * 1.05, floor_value, "MAGeCK RRA floor",
    col = colour_floor, cex = 0.68, pos = 3, adj = 1, xpd = NA
  )
  abline(a = 0, b = 1, col = "#999999", lty = 3)
  points(
    barcs_neglog[!hit], mageck_neglog[!hit], pch = 16, cex = 0.9,
    col = adjustcolor(colour_null, alpha.f = 0.55)
  )
  points(
    barcs_neglog[hit], mageck_neglog[hit], pch = 16, cex = 1.1,
    col = adjustcolor(colour_hit, alpha.f = 0.85)
  )
  to_label <- merged$gene %in% label_genes
  text(
    barcs_neglog[to_label], mageck_neglog[to_label],
    labels = merged$gene[to_label], pos = 4, offset = 0.35, cex = 0.68,
    col = "#333333", xpd = NA
  )
  legend(
    "topleft",
    legend = c(
      sprintf("Spearman rho = %.2f", spearman),
      sprintf("%d/%d MAGeCK genes at the floor", at_floor, nrow(merged)),
      "Called by BARCS (FDR < 0.10)"
    ),
    pch = c(NA, NA, 16), col = c(NA, NA, colour_hit),
    bty = "n", cex = 0.74, text.col = "#333333"
  )
}

draw_figure <- function() {
  layout(matrix(seq_len(2), nrow = 1))
  headtohead_panel(
    -log10(merged$p_value), -log10(merged$rra_p),
    expression(BARCS ~ -log[10](italic(p))),
    expression(MAGeCK ~ RRA ~ -log[10](italic(p))),
    "a", "Significance: -log10 p", -log10(5e-6)
  )
  headtohead_panel(
    -log10(merged$fdr), -log10(merged$rra_fdr),
    expression(BARCS ~ -log[10](FDR)),
    expression(MAGeCK ~ RRA ~ -log[10](FDR)),
    "b", "Significance: -log10 FDR", -log10(min(merged$rra_fdr))
  )
}

dir.create("figures", showWarnings = FALSE)
figure_path <- file.path("figures", "liu_tcell_significance_headtohead.png")
png(figure_path, width = 2000, height = 1050, res = 200)
draw_figure()
invisible(dev.off())

message("Wrote ", figure_path)
message(sprintf(
  "Spearman -log10 p = %.2f, -log10 FDR = %.2f (low values are a MAGeCK RRA tie artefact, not method disagreement)",
  cor(-log10(merged$p_value), -log10(merged$rra_p), method = "spearman"),
  cor(-log10(merged$fdr), -log10(merged$rra_fdr), method = "spearman")
))
