#!/usr/bin/env Rscript

# BARCS versus the published MAGeCK pipeline on the Schmidt et al. primary
# human T cell CRISPRa/CRISPRi screens (Science 2022, doi:10.1126/science.abj4008).
#
# Both methods are run here on the identical GEO raw counts (GSE174255):
#   examples/schmidt_tcell_barcs.R    one beta-binomial model per guide over all
#                                     twelve libraries of its own library set
#   examples/schmidt_tcell_mageck.R   the paper's recipe -- normalise to library
#                                     total, merge the two sets, paired RRA
#
# so every difference below is the model, not the input. BARCS reports a
# logit-scale coefficient and MAGeCK a median log2 fold change; the two are
# compared in sign and rank, never in units.
#
# Writes results/schmidt_tcell/comparison/:
#   <screen>_gene_comparison.csv.gz   both methods' gene-level results, joined
#   <screen>_hits.csv                 the same, cut to genes either method calls
#   method_concordance.csv            per-screen agreement, null scales, hit counts
#   positive_control_panel.csv        rank of each reference gene under each method
#   tcr_pathway_precision.csv         KEGG TCR-signalling enrichment among top-N
# and figures/schmidt_tcell_method_comparison.{png,pdf}. The figure's calibration
# panels are drawn when examples/schmidt_tcell_null_calibration.R has been run.

options(stringsAsFactors = FALSE)

output_dir <- file.path("results", "schmidt_tcell", "output")
mageck_dir <- file.path("results", "schmidt_tcell", "mageck")
comparison_dir <- file.path("results", "schmidt_tcell", "comparison")
figure_dir <- "figures"
dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path("R", "method_palette.R"))
source(file.path("R", "schmidt_tcell.R"))
colour_barcs <- barcs_method_colours[["BARCS"]]
colour_mageck <- barcs_method_colours[["MAGeCK"]]
colour_null <- "#8C8C8C"

screens <- schmidt_screens

# The genes the paper itself names as controls or established regulators. The
# panel is fixed from the paper's text before either method is looked at; it is
# the reference both are scored against.
#
# Only each gene's biological role is asserted. The expected sort-bin direction
# follows from the perturbation: activating a positive regulator of the sorted
# cytokine enriches the guide in the high bin, silencing it depletes the guide,
# and a negative regulator does the reverse.
reference_roles <- rbind(
  # Screen positive controls: the sorted cytokine's own gene (Fig. 1B).
  data.frame(gene = "IL2",  screen = c("CRISPRa_IL2", "CRISPRi_IL2"),
             role = "positive control", promotes = TRUE),
  data.frame(gene = "IFNG", screen = c("CRISPRa_IFNG", "CRISPRi_IFNG"),
             role = "positive control", promotes = TRUE),
  # "VAV1, CD28, LCP2 (encoding SLP-76), and LAT reinforced T cell activation
  # and were enriched in both cytokine-high bins", with ZAP70 named as a core
  # positive regulator across all four screens.
  expand.grid(gene = c("VAV1", "CD28", "LCP2", "LAT", "ZAP70"), screen = screens,
              role = "proximal TCR signalling", promotes = TRUE,
              stringsAsFactors = FALSE),
  # "the negative TCR signalling regulators MAP4K1 and SLA2 were depleted in
  # these bins", with CBLB named as the core negative regulator.
  expand.grid(gene = c("MAP4K1", "SLA2", "CBLB"), screen = screens,
              role = "negative regulator", promotes = FALSE,
              stringsAsFactors = FALSE),
  # "CRISPRi detected a circuit of T cell stimulation signalling through MALT1,
  # BCL10, TRAF6, and TAK1 (encoded by MAP3K7) to the inhibitor of the NF-kB
  # complex (encoded by CHUK, IKBKB, and IKBKG) that promotes IFN-g production"
  # -- asserted for CRISPRi and IFN-g only.
  expand.grid(gene = c("MALT1", "BCL10", "TRAF6", "MAP3K7", "CHUK", "IKBKB", "IKBKG"),
              screen = "CRISPRi_IFNG", role = "NF-kB circuit", promotes = TRUE,
              stringsAsFactors = FALSE)
)
reference_panel <- transform(
  reference_roles,
  gene = as.character(gene),
  screen = as.character(screen),
  # CRISPRa raises the target, CRISPRi lowers it, so the two modalities expect
  # opposite bin shifts from the same biological role.
  direction = ifelse(promotes, 1, -1) *
    ifelse(startsWith(as.character(screen), "CRISPRa"), 1, -1)
)

read_barcs <- function(screen) {
  path <- file.path(output_dir, paste0(screen, "_GENES.csv.gz"))
  if (!file.exists(path)) {
    stop("Missing ", path, ". Run `Rscript examples/schmidt_tcell_barcs.R` first.",
         call. = FALSE)
  }
  genes <- read.csv(path)
  # The nontargeting pseudo-genes are a diagnostic, not screen results.
  genes <- genes[!genes$control_gene, ]
  # `bb_gene_consistency()` scales by max(all-gene MAD, control-derived scale),
  # so whichever of the two happens to be larger sets the scale for the whole
  # screen. This variant -- same fits, control-derived scale alone -- isolates
  # what that rule costs. examples/schmidt_tcell_null_calibration.R checks on
  # held-out controls whether the extra calls survive.
  control_statistic <- (genes$raw_statistic - genes$null_center) /
    genes$control_scale
  data.frame(
    gene = genes$gene,
    n_guides = genes$n_guides,
    barcs_estimate = genes$estimate,
    barcs_p = genes$p_value,
    barcs_fdr = genes$fdr,
    barcs_agreement = genes$guide_direction_agreement,
    barcs_null_scale = genes$null_scale,
    barcs_control_scale = genes$control_scale,
    barcs_control_scale_fdr = p.adjust(
      2 * pnorm(-abs(control_statistic)), method = "BH"
    )
  )
}

read_mageck <- function(screen) {
  path <- file.path(mageck_dir, screen, paste0(screen, ".gene_summary.txt"))
  if (!file.exists(path)) {
    stop("Missing ", path, ". Run `Rscript examples/schmidt_tcell_mageck.R` first.",
         call. = FALSE)
  }
  genes <- read.delim(path, check.names = FALSE)
  # RRA is one-sided in each direction. The published hit rule reads a single
  # FDR per gene, which is the smaller of the two tails.
  #
  # Ranking uses the RRA rho score rather than the permutation p-value: the
  # p-value saturates -- 124 genes tie at its minimum in CRISPRa IL-2 -- and
  # ranking on a saturated statistic would understate MAGeCK's ordering.
  data.frame(
    gene = genes$id,
    mageck_lfc = genes[["pos|lfc"]],
    mageck_score = pmin(genes[["pos|score"]], genes[["neg|score"]]),
    mageck_p = pmin(genes[["pos|p-value"]], genes[["neg|p-value"]]),
    mageck_fdr = pmin(genes[["pos|fdr"]], genes[["neg|fdr"]]),
    mageck_enriched = genes[["pos|p-value"]] <= genes[["neg|p-value"]]
  )
}

# `schmidt_mageck_hits()` applied to the joined table rather than to a raw
# gene_summary: the paper's rule, median |log2FC| > 0.5 and FDR < 0.05.
mageck_hits <- function(x) abs(x$mageck_lfc) > 0.5 & x$mageck_fdr < 0.05

comparisons <- list()
for (screen in screens) {
  joined <- merge(read_barcs(screen), read_mageck(screen), by = "gene")
  joined <- joined[joined$gene != "NO-TARGET", ]
  joined$barcs_rank <- rank(joined$barcs_p, ties.method = "min")
  joined$mageck_rank <- rank(joined$mageck_score, ties.method = "min")
  comparisons[[screen]] <- joined
  write.csv(
    joined, gzfile(file.path(comparison_dir, paste0(screen, "_gene_comparison.csv.gz"))),
    row.names = FALSE
  )
  # The full 18,800-gene join is regenerable output. What is worth carrying in
  # the repository is every gene either method calls, plus the reference panel:
  # a few hundred rows that support every number quoted in the write-up.
  keep <- joined$barcs_fdr < 0.05 | mageck_hits(joined) |
    joined$gene %in% reference_panel$gene[reference_panel$screen == screen]
  hits <- joined[which(keep), ]
  hits <- hits[order(pmin(hits$barcs_rank, hits$mageck_rank)), ]
  write.csv(hits, file.path(comparison_dir, paste0(screen, "_hits.csv")),
            row.names = FALSE)
}

# ---- per-screen concordance and hit counts -------------------------------

concordance <- do.call(rbind, lapply(screens, function(screen) {
  x <- comparisons[[screen]]
  barcs_hit <- x$barcs_fdr < 0.05
  mageck_hit <- mageck_hits(x)
  data.frame(
    screen = screen,
    genes = nrow(x),
    spearman_effect = cor(x$barcs_estimate, x$mageck_lfc, method = "spearman",
                          use = "complete.obs"),
    sign_agreement = mean(sign(x$barcs_estimate) == sign(x$mageck_lfc), na.rm = TRUE),
    barcs_null_scale = x$barcs_null_scale[1L],
    barcs_control_scale = x$barcs_control_scale[1L],
    # Hit counts under three rules. The published MAGeCK rule is a conjunction
    # -- FDR < 0.05 *and* median |log2FC| > 0.5 -- and the effect floor is not
    # cosmetic: it removes 38-54% of the FDR-significant genes and halves the
    # held-out null false positives. Comparing BARCS at a bare FDR against that
    # conjunction compares two different rule shapes, so the bare-FDR MAGeCK
    # count is reported alongside it as the like-for-like statistical rule.
    barcs_hits = sum(barcs_hit, na.rm = TRUE),
    mageck_hits = sum(mageck_hit, na.rm = TRUE),
    mageck_hits_fdr_only = sum(x$mageck_fdr < 0.05, na.rm = TRUE),
    barcs_control_scale_hits = sum(x$barcs_control_scale_fdr < 0.05, na.rm = TRUE),
    shared_hits = sum(barcs_hit & mageck_hit, na.rm = TRUE),
    spearman_top1000 = local({
      top <- x[x$barcs_rank <= 1000 | x$mageck_rank <= 1000, ]
      cor(top$barcs_estimate, top$mageck_lfc, method = "spearman", use = "complete.obs")
    })
  )
}))
write.csv(concordance, file.path(comparison_dir, "method_concordance.csv"), row.names = FALSE)

# ---- reference panel: where each method ranks the known answers -----------

panel <- do.call(rbind, lapply(seq_len(nrow(reference_panel)), function(i) {
  entry <- reference_panel[i, ]
  x <- comparisons[[entry$screen]]
  row <- x[match(entry$gene, x$gene), ]
  if (is.na(row$gene)) return(NULL)
  data.frame(
    screen = entry$screen, gene = entry$gene, role = entry$role,
    expected_direction = entry$direction,
    barcs_rank = row$barcs_rank, barcs_estimate = row$barcs_estimate,
    barcs_fdr = row$barcs_fdr,
    barcs_called = row$barcs_fdr < 0.05 &&
      sign(row$barcs_estimate) == entry$direction,
    mageck_rank = row$mageck_rank, mageck_lfc = row$mageck_lfc,
    mageck_fdr = row$mageck_fdr,
    mageck_called = mageck_hits(row) && sign(row$mageck_lfc) == entry$direction
  )
}))
write.csv(panel, file.path(comparison_dir, "positive_control_panel.csv"), row.names = FALSE)

# ---- precision proxy: KEGG TCR signalling among each method's top genes ----
#
# Ranking genes is the part of a screen analysis that matters, and the two
# methods return different orderings of the same data. A curated pathway the
# screens are known to act through gives a label-free precision proxy: at a
# matched list length, the method whose top genes are more enriched for KEGG T
# cell receptor signalling is ordering the screen better. This is a proxy, not
# a truth set -- most real hits lie outside any one pathway.
tcr_precision <- NULL
if (requireNamespace("msigdbr", quietly = TRUE)) {
  # msigdbr renamed its arguments and the KEGG subcollection at version 10.
  kegg <- try(
    msigdbr::msigdbr(species = "Homo sapiens", collection = "C2",
                     subcollection = "CP:KEGG_LEGACY"),
    silent = TRUE
  )
  if (inherits(kegg, "try-error")) {
    kegg <- try(
      msigdbr::msigdbr(species = "Homo sapiens", category = "C2",
                       subcategory = "CP:KEGG"),
      silent = TRUE
    )
  }
  if (!inherits(kegg, "try-error")) {
    is_tcr <- kegg$gs_name == "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY"
    tcr_genes <- unique(kegg$gene_symbol[is_tcr])
    tcr_precision <- do.call(rbind, lapply(screens, function(screen) {
      x <- comparisons[[screen]]
      background <- mean(x$gene %in% tcr_genes)
      do.call(rbind, lapply(c(100L, 250L, 500L), function(n) {
        data.frame(
          screen = screen, top_n = n,
          background_rate = background,
          barcs_tcr = sum(x$gene[x$barcs_rank <= n] %in% tcr_genes),
          mageck_tcr = sum(x$gene[x$mageck_rank <= n] %in% tcr_genes)
        )
      }))
    }))
    tcr_precision$barcs_fold <- (tcr_precision$barcs_tcr / tcr_precision$top_n) /
      tcr_precision$background_rate
    tcr_precision$mageck_fold <- (tcr_precision$mageck_tcr / tcr_precision$top_n) /
      tcr_precision$background_rate
    write.csv(tcr_precision, file.path(comparison_dir, "tcr_pathway_precision.csv"),
              row.names = FALSE)
  }
}

print(concordance, row.names = FALSE)
cat("\nReference panel, genes called by each method:\n")
print(with(panel, table(role, BARCS = barcs_called, MAGeCK = mageck_called)))
cat("\nPanel called under each method's own rule: BARCS", sum(panel$barcs_called),
    "/", nrow(panel), " MAGeCK", sum(panel$mageck_called), "/", nrow(panel), "\n")
cat("Panel recovered at matched list length:\n")
for (n in c(50L, 100L, 250L, 500L)) {
  cat(sprintf("  top %4d: BARCS %2d/%d  MAGeCK %2d/%d\n", n,
              sum(panel$barcs_rank <= n), nrow(panel),
              sum(panel$mageck_rank <= n), nrow(panel)))
}
if (!is.null(tcr_precision)) {
  cat("\nKEGG TCR-signalling genes among each method's top N:\n")
  print(tcr_precision, row.names = FALSE)
}

# ---- figure ---------------------------------------------------------------

# The held-out null calibration is produced by
# examples/schmidt_tcell_null_calibration.R, which needs MAGeCK. Its panels are
# drawn when that file is present.
null_path <- file.path(comparison_dir, "null_calibration.csv")
null_calibration <- if (file.exists(null_path)) read.csv(null_path) else NULL

method_colours <- c(
  BARCS = colour_barcs,
  `BARCS unmoderated` = adjustcolor(colour_barcs, alpha.f = 0.45),
  MAGeCK = colour_mageck
)

draw <- function() {
  op <- par(mfrow = c(2, 4), mar = c(4.2, 4.2, 2.6, 1.0), mgp = c(2.5, 0.8, 0),
            cex.axis = 0.8, cex.lab = 0.9, cex.main = 1.0)
  on.exit(par(op), add = TRUE)

  # (a-d) effect concordance, one panel per screen
  for (screen in screens) {
    x <- comparisons[[screen]]
    hit <- x$barcs_fdr < 0.05 | mageck_hits(x)
    plot(x$mageck_lfc, x$barcs_estimate, pch = 16, cex = 0.28,
         col = ifelse(hit, colour_barcs, adjustcolor(colour_null, alpha.f = 0.25)),
         xlab = "MAGeCK median log2(high/low)", ylab = "BARCS effect (logit)",
         main = sub("_", " ", screen))
    abline(h = 0, v = 0, col = colour_null, lty = 3)
    legend("topleft", bty = "n", cex = 0.75,
           legend = sprintf("Spearman %.2f", cor(x$barcs_estimate, x$mageck_lfc,
                                                 method = "spearman", use = "complete.obs")))
  }

  # (e) where each method ranks the genes the paper names
  ordered <- panel[order(panel$screen, panel$barcs_rank), ]
  par(mar = c(4.2, 6.4, 2.6, 0.6))
  plot(NA, xlim = range(c(ordered$barcs_rank, ordered$mageck_rank)),
       ylim = c(nrow(ordered) + 0.5, 0.5), log = "x", yaxt = "n",
       xlab = "rank in screen (log)", ylab = "", main = "Genes the paper names")
  rows <- seq_len(nrow(ordered))
  segments(ordered$barcs_rank, rows, ordered$mageck_rank, rows, col = colour_null)
  points(ordered$barcs_rank, rows, pch = 16, cex = 0.55, col = colour_barcs)
  points(ordered$mageck_rank, rows, pch = 17, cex = 0.55, col = colour_mageck)
  axis(2, at = rows, labels = ordered$gene, las = 1, cex.axis = 0.4, tick = FALSE,
       line = -0.6)
  # One rule per screen block, so the four screens read apart.
  breaks <- which(diff(match(ordered$screen, screens)) != 0) + 0.5
  abline(h = breaks, col = colour_null, lty = 3)
  # The x axis is logarithmic, so `par("usr")` is in log10 units.
  text(10^par("usr")[2L], c(0.5, breaks) + 1.6, adj = c(1, 0.5), cex = 0.5,
       col = colour_null, labels = sub("_", " ", screens))
  legend("bottomleft", bty = "n", cex = 0.7, pch = c(16, 17),
         col = c(colour_barcs, colour_mageck), legend = c("BARCS", "MAGeCK"))
  par(mar = c(4.2, 4.2, 2.6, 1.0))

  # (f) hit counts. The published MAGeCK rule is a conjunction -- FDR < 0.05 and
  # median |log2FC| > 0.5 -- so it is shown next to MAGeCK at a bare FDR, which
  # is the rule shape BARCS is read at.
  counts <- rbind(concordance$barcs_hits, concordance$mageck_hits_fdr_only,
                  concordance$mageck_hits)
  colnames(counts) <- sub("_", "\n", screens)
  barplot(counts, beside = TRUE,
          col = c(colour_barcs, adjustcolor(colour_mageck, alpha.f = 0.45),
                  colour_mageck),
          border = NA, ylab = "genes called", main = "Hits per screen",
          cex.names = 0.7,
          legend.text = c("BARCS, FDR<0.05", "MAGeCK, FDR<0.05",
                          "MAGeCK, published rule"),
          args.legend = list(bty = "n", cex = 0.6, x = "topleft"))

  if (is.null(null_calibration)) {
    plot.new()
    plot.new()
    return(invisible(NULL))
  }

  # (g) held-out nontargeting pseudo-genes: observed rate below nominal 0.05,
  # against each method's matched expectation. BARCS reports one two-sided
  # p-value, so 5%; the published MAGeCK rule reads the smaller of two
  # near-complementary one-sided RRA tails, so 10%.
  ratio <- with(null_calibration, tapply(null_p05_ratio, list(method, screen), mean))
  ratio <- ratio[names(method_colours), screens, drop = FALSE]
  colnames(ratio) <- sub("_", "\n", screens)
  barplot(ratio, beside = TRUE, col = method_colours, border = NA,
          ylab = "observed / expected null rate", cex.names = 0.7,
          main = "Held-out null calibration")
  abline(h = 1, col = colour_null, lty = 2)

  # (h) what the null scale costs in stability: one point per held-out fit.
  spread <- split(null_calibration$genes_called, null_calibration$method)
  spread <- spread[names(method_colours)]
  plot(NA, xlim = c(0.5, length(spread) + 0.5),
       ylim = range(null_calibration$genes_called), xaxt = "n",
       xlab = "", ylab = "genes called", main = "Stability over held-out fits")
  axis(1, at = seq_along(spread),
       labels = c("BARCS", "un-\nmoderated", "MAGeCK"),
       cex.axis = 0.7, padj = 0.5)
  for (i in seq_along(spread)) {
    points(jitter(rep(i, length(spread[[i]])), amount = 0.12), spread[[i]],
           pch = 16, cex = 0.8, col = method_colours[[i]])
    segments(i - 0.28, median(spread[[i]]), i + 0.28, median(spread[[i]]),
             col = method_colours[[i]], lwd = 2)
  }
}

png(file.path(figure_dir, "schmidt_tcell_method_comparison.png"),
    width = 2800, height = 1400, res = 200)
draw()
invisible(dev.off())
pdf(file.path(figure_dir, "schmidt_tcell_method_comparison.pdf"),
    width = 14.0, height = 7.0)
draw()
invisible(dev.off())
message("wrote figures/schmidt_tcell_method_comparison.{png,pdf}")
