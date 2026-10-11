#!/usr/bin/env Rscript

# BARCS versus the published MAGeCK analysis of the genome-wide Liu et al.
# in vivo T cell IFNγ screen (Nature 2026, doi:10.1038/s41586-026-10906-9,
# Supplementary Table 2, the screen behind main Fig. 4).
#
# The paper summarises the screen with a one-sided MAGeCK RRA test and a
# per-gene log2 fold change between the IFNγ-high and IFNγ-low sort gates.
# BARCS reports a logit-scale beta-binomial coefficient with a calibrated
# p-value. The two are comparable in sign and rank, not in units.
#
# Panels:
#   a  BARCS effect vs MAGeCK gene log2FC across all ~19k genes (concordance;
#      cf. main Fig. 4b). Genes called by BARCS at FDR < 0.10 are highlighted.
#   b  BARCS genome-wide volcano (effect vs -log10 p), the BARCS counterpart of
#      the MAGeCK RRA volcano in Fig. 4b.
#
# Requires results/liu_genomewide/output/ (run examples/liu_genomewide_barcs.R)
# and the workbook under data/raw/liu_tcell/ (run scripts/prepare_liu_genomewide.R).

options(stringsAsFactors = FALSE)

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("The `readxl` package is required.")
}

output_dir <- file.path("results", "liu_genomewide", "output")
table2_path <- file.path("data", "raw", "liu_tcell", "Supplementary_Table_2.xlsx")
genes_path <- file.path(output_dir, "gate-donor_GENES.csv")
if (!file.exists(genes_path)) {
  stop("Missing ", genes_path, ". Run `Rscript examples/liu_genomewide_barcs.R` first.",
       call. = FALSE)
}
if (!file.exists(table2_path)) {
  stop("Missing ", table2_path, ". Run `Rscript scripts/prepare_liu_genomewide.R` first.",
       call. = FALSE)
}

barcs <- read.csv(genes_path)
barcs <- barcs[barcs$gene != "NTCTRL", ]

published <- as.data.frame(
  suppressMessages(readxl::read_excel(table2_path, sheet = "IFNghi vs. IFNglo_gene_summary")),
  check.names = FALSE
)
# MAGeCK reports one log2FC per gene in both blocks; the RRA p-value is the
# significant tail.
published <- data.frame(
  gene = published$id,
  log2fc = published[["pos|lfc"]],
  rra_fdr = pmin(published[["pos|fdr"]], published[["neg|fdr"]])
)

merged <- merge(
  barcs[, c("gene", "estimate", "p_value", "fdr")], published, by = "gene"
)

hit_fdr <- 0.10
colour_hit <- "#0072B2"
colour_null <- "#8C8C8C"
colour_enriched <- "#D55E00"

# Label only the canonical IFNγ controls plus the few strongest hits: on a
# 19k-gene scatter, dense labelling collides into noise rather than informing.
control_genes <- c("IFNG", "TNF", "STAT1")
top_barcs <- head(merged$gene[order(merged$fdr, -abs(merged$estimate))], 5)
label_genes <- unique(c(control_genes, top_barcs))

concordance_panel <- function() {
  hit <- merged$fdr < hit_fdr
  spearman <- cor(merged$estimate, merged$log2fc, method = "spearman")
  agreement <- mean(sign(merged$estimate) == sign(merged$log2fc))

  par(mar = c(4.2, 4.4, 2.8, 1.2), cex.axis = 0.92, cex.lab = 0.98)
  plot(
    merged$estimate, merged$log2fc, type = "n",
    xlab = "BARCS effect (logit-scale coefficient)",
    ylab = "MAGeCK gene log2 fold change",
    main = "a  Genome-wide concordance (cf. Fig. 4b)", bty = "l", cex.main = 1.02
  )
  abline(h = 0, col = "#666666", lty = 3)
  abline(v = 0, col = "#666666", lty = 3)
  points(
    merged$estimate[!hit], merged$log2fc[!hit], pch = 16, cex = 0.35,
    col = adjustcolor(colour_null, alpha.f = 0.30)
  )
  points(
    merged$estimate[hit], merged$log2fc[hit], pch = 16, cex = 0.55,
    col = adjustcolor(colour_hit, alpha.f = 0.60)
  )
  to_label <- merged$gene %in% label_genes
  text(
    merged$estimate[to_label], merged$log2fc[to_label],
    labels = merged$gene[to_label], pos = 3, offset = 0.3, cex = 0.62,
    col = "#333333", xpd = NA
  )
  legend(
    "topleft",
    legend = c(
      sprintf("Spearman rho = %.2f", spearman),
      sprintf("sign agreement = %.0f%%", 100 * agreement),
      sprintf("%d genes; blue = BARCS FDR < 0.10", nrow(merged))
    ),
    pch = c(NA, NA, 16), col = c(NA, NA, colour_hit),
    bty = "n", cex = 0.72, text.col = "#333333"
  )
}

volcano_panel <- function() {
  effect <- merged$estimate
  y <- -log10(merged$p_value)
  hit <- merged$fdr < hit_fdr
  base_colour <- ifelse(!hit, colour_null,
                        ifelse(effect > 0, colour_enriched, colour_hit))
  alpha <- ifelse(hit, 0.6, 0.3)
  point_colour <- mapply(function(c, a) adjustcolor(c, alpha.f = a),
                         base_colour, alpha)

  par(mar = c(4.2, 4.4, 2.8, 1.2), cex.axis = 0.92, cex.lab = 0.98)
  plot(
    effect, y, type = "n",
    ylim = c(0, max(y) * 1.15),
    xlab = "BARCS effect (logit-scale coefficient)",
    ylab = expression(-log[10](italic(p))),
    main = "b  BARCS genome-wide volcano (cf. RRA Fig. 4b)", bty = "l",
    cex.main = 1.02
  )
  abline(v = 0, col = "#666666", lty = 3)
  abline(h = -log10(0.05), col = "#3B5BDB", lty = 2)
  points(effect, y, pch = 16, cex = ifelse(hit, 0.55, 0.35), col = point_colour)
  to_label <- merged$gene %in% label_genes & hit
  text(
    effect[to_label], y[to_label], labels = merged$gene[to_label],
    pos = 3, offset = 0.3, cex = 0.62, col = "#333333", xpd = NA
  )
  legend(
    "topright",
    legend = c("Enriched (FDR < 0.10)", "Depleted (FDR < 0.10)", "Not called"),
    pch = 16, cex = 0.72, bty = "n", text.col = "#333333",
    col = c(colour_enriched, colour_hit, colour_null)
  )
}

draw_figure <- function() {
  layout(matrix(seq_len(2), nrow = 1))
  concordance_panel()
  volcano_panel()
}

dir.create("figures", showWarnings = FALSE)
figure_stem <- file.path("figures", "liu_genomewide_publication_comparison")
png(paste0(figure_stem, ".png"), width = 2000, height = 1050, res = 200)
draw_figure()
invisible(dev.off())
pdf(paste0(figure_stem, ".pdf"), width = 10, height = 5.25, pointsize = 11,
    useDingbats = FALSE)
draw_figure()
invisible(dev.off())

message("Wrote ", figure_stem, ".png and ", figure_stem, ".pdf")
message(sprintf(
  "Spearman(effect, MAGeCK log2FC) = %.2f over %d genes",
  cor(merged$estimate, merged$log2fc, method = "spearman"), nrow(merged)
))
