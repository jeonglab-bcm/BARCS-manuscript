#!/usr/bin/env Rscript

# Liang HAP1 processed-count example modeled on a specificity-threshold curve
# plus side-by-side volcano plots. BARCS, official MAGeCK-MLE, edgeR-QL,
# DESeq2, and limma-voom use the same rounded normalized/ComBat-corrected
# pseudo-count matrix and the same replicate/day design.
#
# Every other real-data benchmark here measures whether a method finds the
# genes it should. This one measures whether it stays quiet when it should,
# which is the failure mode that costs a screen the most: a false hit is
# followed up at bench cost, and nothing in the screen itself reveals the
# error.
#
# HAP1 supplies null controls, so the fraction of them called at a given
# threshold is a direct false-positive rate rather than an estimate. Sweeping
# the threshold turns that into a curve, and the curve is the honest way to
# compare methods whose p-values are on different scales -- comparing them at
# one arbitrary cutoff would mostly compare their calibration.
#
# The volcano panels beside it exist so the curve can be interpreted. Two
# methods can share a false-positive rate and get there differently, one by
# shrinking everything toward the middle and one by being genuinely
# discriminating, and only the volcano shows which.
#
# All five methods are handed the identical matrix and design. This is
# deliberate and it is a limitation: the input is Liang's already normalized,
# ComBat-corrected, rounded pseudo-counts, not raw counts, so the count models
# are being run outside the input they were designed for. It equalizes
# preprocessing at the cost of flattering none of them.
#
#     Rscript examples/liang_cas13_benchmark.R    # writes the scores read here
#     Rscript examples/liang_hap1_specificity_volcano.R

options(stringsAsFactors = FALSE)
source(file.path("R", "method_palette.R"))
analysis_protocol <- "liang-hap1-specificity-volcano-v2"
score_path <- file.path(
  "results", "liang_cas13", "all_gene_scores.csv.gz"
)
if (!file.exists(score_path)) {
  stop(
    "Run `examples/liang_cas13_benchmark.R` before this example.",
    call. = FALSE
  )
}

methods <- c(
  "BARCS", "MAGeCK-MLE", "edgeR-QL", "DESeq2", "limma-voom"
)
thresholds <- c(0.001, 0.005, 0.01, 0.05, 0.10, 0.20)
scores <- read.csv(gzfile(score_path))
hap1 <- scores[
  scores$cell_line == "HAP1" &
    scores$method %in% methods &
    is.finite(scores$effect) &
    is.finite(scores$p_value) &
    is.finite(scores$fdr),
  ,
  drop = FALSE
]
if (!identical(sort(unique(hap1$method)), sort(methods))) {
  stop(
    "The HAP1 BARCS, MAGeCK-MLE, edgeR-QL, DESeq2, and ",
    "limma-voom results are required."
  )
}
if (!all(
  hap1$input_scale ==
    "rounded_normalized_ComBat_pseudocount"
)) {
  stop("Unexpected Liang input scale.")
}

specificity_rows <- list()
row_index <- 0L
for (method in methods) {
  method_data <- hap1[hap1$method == method, , drop = FALSE]
  null_data <- method_data[
    !is.na(method_data$truth) & method_data$truth == 0,
    ,
    drop = FALSE
  ]
  if (nrow(null_data) < 1500L) {
    stop("Too few HAP1 non-expressed lncRNA null controls for ", method)
  }
  for (threshold in thresholds) {
    row_index <- row_index + 1L
    false_positive <- sum(null_data$p_value < threshold)
    specificity_rows[[row_index]] <- data.frame(
      analysis_protocol = analysis_protocol,
      cell_line = "HAP1",
      input_scale = unique(method_data$input_scale),
      method = method,
      p_value_threshold = threshold,
      null_controls = nrow(null_data),
      false_positive_nulls = false_positive,
      specificity = 1 - false_positive / nrow(null_data)
    )
  }
}
specificity <- do.call(rbind, specificity_rows)
write.csv(
  specificity,
  file.path(
    "data", "derived",
    "liang_hap1_specificity_by_p_threshold.csv"
  ),
  row.names = FALSE
)

volcano_summary <- do.call(rbind, lapply(methods, function(method) {
  method_data <- hap1[hap1$method == method, , drop = FALSE]
  null_data <- method_data[
    !is.na(method_data$truth) & method_data$truth == 0,
    ,
    drop = FALSE
  ]
  essential_data <- method_data[
    !is.na(method_data$truth) & method_data$truth == 1,
    ,
    drop = FALSE
  ]
  data.frame(
    analysis_protocol = analysis_protocol,
    cell_line = "HAP1",
    input_scale = unique(method_data$input_scale),
    method = method,
    genes = nrow(method_data),
    discoveries_fdr_0_10 = sum(method_data$fdr < 0.10),
    depleted_discoveries_fdr_0_10 =
      sum(method_data$fdr < 0.10 & method_data$effect < 0),
    essential_genes = nrow(essential_data),
    essential_depleted_recovered_fdr_0_10 = sum(
      essential_data$fdr < 0.10 & essential_data$effect < 0
    ),
    null_controls = nrow(null_data),
    null_p_below_0_05 = mean(null_data$p_value < 0.05)
  )
}))
write.csv(
  volcano_summary,
  file.path(
    "data", "derived", "liang_hap1_volcano_summary.csv"
  ),
  row.names = FALSE
)

method_colours <- c(
  BARCS = barcs_method_colours[["BARCS"]],
  `MAGeCK-MLE` = barcs_method_colours[["MAGeCK"]],
  `edgeR-QL` = barcs_method_colours[["edgeR"]],
  DESeq2 = barcs_method_colours[["DESeq2"]],
  `limma-voom` = barcs_method_colours[["limma-voom"]]
)
method_pch <- c(
  BARCS = 16, `MAGeCK-MLE` = 17, `edgeR-QL` = 15,
  DESeq2 = 18, `limma-voom` = 8
)

pdf(
  file.path("figures", "liang_hap1_specificity_volcano.pdf"),
  width = 12,
  height = 7.6,
  useDingbats = FALSE
)
layout(matrix(1:6, nrow = 2, byrow = TRUE))

par(mar = c(5.2, 4.5, 3, 1.2))
plot(
  NA,
  xlim = range(seq_along(thresholds)),
  ylim = c(
    min(specificity$specificity) - 0.03,
    min(1, max(specificity$specificity) + 0.015)
  ),
  xaxt = "n",
  xlab = "Nominal p-value threshold",
  ylab = "Specificity among null lncRNAs",
  main = "A  HAP1 null specificity",
  bty = "l"
)
axis(
  1,
  at = seq_along(thresholds),
  labels = format(thresholds, trim = TRUE, scientific = FALSE),
  las = 2,
  cex.axis = 0.85
)
for (method in methods) {
  panel <- specificity[specificity$method == method, ]
  panel <- panel[
    match(thresholds, panel$p_value_threshold),
    ,
    drop = FALSE
  ]
  lines(
    seq_along(thresholds),
    panel$specificity,
    type = "o",
    pch = method_pch[[method]],
    lwd = 2,
    col = method_colours[[method]]
  )
}
legend(
  "bottomleft",
  legend = methods,
  col = unname(method_colours[methods]),
  pch = unname(method_pch[methods]),
  lty = 1,
  lwd = 2,
  bty = "n",
  cex = 0.72
)

common_xlim <- range(hap1$effect)
y_cap <- 60
common_ylim <- c(0, y_cap)
panel_letters <- LETTERS[2:6]
for (method_index in seq_along(methods)) {
  method <- methods[method_index]
  panel <- hap1[hap1$method == method, , drop = FALSE]
  y_value <- -log10(pmax(panel$p_value, .Machine$double.xmin))
  y_clipped <- y_value > y_cap
  depleted_call <- panel$fdr < 0.10 & panel$effect < 0
  point_colours <- ifelse(
    depleted_call,
    adjustcolor("#E64B35", alpha.f = 0.55),
    adjustcolor("#999999", alpha.f = 0.25)
  )
  show_y_axis <- method_index %in% c(1L, 3L)
  par(mar = c(5.2, if (show_y_axis) 4.5 else 1.2, 3, 1))
  plot(
    panel$effect,
    pmin(y_value, y_cap),
    pch = ifelse(y_clipped, 17, 16),
    cex = 0.38,
    col = point_colours,
    xlim = common_xlim,
    ylim = common_ylim,
    xlab = "Gene effect",
    ylab = if (show_y_axis) {
      expression(-log[10](italic(p))~"(cap 60)")
    } else {
      ""
    },
    yaxt = if (show_y_axis) "s" else "n",
    main = paste0(panel_letters[method_index], "  ", method),
    bty = "l"
  )
  abline(h = -log10(0.05), col = "#3B5BDB", lty = 2)
  abline(v = 0, col = "#666666", lty = 3)
}
dev.off()

print(specificity)
print(volcano_summary)
