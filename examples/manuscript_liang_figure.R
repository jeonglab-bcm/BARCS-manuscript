#!/usr/bin/env Rscript

# Manuscript Figure 1: longitudinal Liang method comparison.
#
# Panels A-E show the four replicate-complete cell lines in method-specific
# volcano plots for BARCS, MAGeCK-MLE, edgeR-QL, DESeq2, and limma-voom, on a
# shared x and y scale. Panels F-H show the observed day 0, 7, and 14
# guide-abundance trajectories for representative disagreements within the
# prespecified proxy-null set. Colour means cell line throughout.

options(stringsAsFactors = FALSE)
source(file.path("R", "method_palette.R"))
source(file.path("R", "load_barcs.R"))

result_dir <- file.path("results", "liang_cas13")
score_path <- file.path(result_dir, "all_gene_scores.csv.gz")
if (!file.exists(score_path)) {
  stop(
    "Run `examples/liang_cas13_benchmark.R` before drawing Figure 1.",
    call. = FALSE
  )
}

cell_stems <- c(
  HAP1 = "HAP1",
  HEK293FT = "HEK293FT",
  `MDA-MB-231` = "MDA_MB_231",
  THP1 = "THP1"
)
guide_method_files <- c(
  BARCS = "barcs",
  DESeq2 = "deseq2",
  `edgeR-QL` = "edger_ql",
  `limma-voom` = "limma_voom"
)
case_specification <- data.frame(
  cell_line = c("HAP1", "MDA-MB-231", "HEK293FT"),
  gene = c(
    "Hum_XLOC_013346",
    "Hum_XLOC_029247",
    "Hum_XLOC_000540"
  ),
  guide = c("gL_036790", "gL_017527", "gL_000389"),
  stringsAsFactors = FALSE
)

gene_scores <- read.csv(gzfile(score_path))
guide_truth <- unique(gene_scores[
  gene_scores$method == "BARCS",
  c("cell_line", "gene", "truth")
])

require_result_file <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Missing ", path, ". Run `examples/liang_cas13_benchmark.R` before ",
      "drawing Figure 1.",
      call. = FALSE
    )
  }
  path
}

load_guide_results <- function(cell_line, method) {
  cell_stem <- cell_stems[[cell_line]]
  path <- require_result_file(file.path(
    result_dir,
    paste0(
      cell_stem, "_longitudinal_",
      guide_method_files[[method]], "_guide.csv.gz"
    )
  ))
  result <- read.csv(gzfile(path))
  is_control <- result$gene == "non-targeting" |
    startsWith(result$gene, "NT_")
  signed_z <- sign(result$estimate) * qnorm(
    pmax(result$p_value / 2, .Machine$double.xmin),
    lower.tail = FALSE
  )
  scale <- max(
    1,
    unname(quantile(abs(signed_z[is_control]), 0.95, type = 8)) /
      qnorm(0.975)
  )
  result$p_value <- 2 * pnorm(-abs(signed_z / scale))
  result
}

load_case_results <- function(cell_line, gene, guide) {
  truth_row <- guide_truth[
    guide_truth$cell_line == cell_line &
      guide_truth$gene == gene,
    ,
    drop = FALSE
  ]
  if (nrow(truth_row) != 1L || truth_row$truth != 0) {
    stop("Single-guide example is not a prespecified null: ", cell_line, "/", gene)
  }

  result <- lapply(names(guide_method_files), function(method) {
    table <- load_guide_results(cell_line, method)
    row <- table[
      table$gene == gene & table$guide == guide,
      ,
      drop = FALSE
    ]
    if (nrow(row) != 1L) {
      stop("Could not resolve one guide row for ", cell_line, "/", guide)
    }
    data.frame(
      cell_line = cell_line,
      gene = gene,
      guide = guide,
      method = method,
      effect = row$estimate,
      p_value = row$p_value,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, result)
}

case_results <- do.call(rbind, lapply(
  seq_len(nrow(case_specification)),
  function(index) {
    load_case_results(
      case_specification$cell_line[index],
      case_specification$gene[index],
      case_specification$guide[index]
    )
  }
))
if (any(
  case_results$p_value[case_results$method == "BARCS"] < 0.20
) || any(
  case_results$p_value[case_results$method != "BARCS"] >= 0.05
)) {
  stop("Single-guide examples no longer satisfy the stated selection rule.")
}
write.csv(
  case_results,
  file.path(
    "data", "derived", "liang_single_guide_null_examples.csv"
  ),
  row.names = FALSE
)

guide_disagreement <- do.call(rbind, lapply(
  names(cell_stems),
  function(cell_line) {
    method_tables <- lapply(names(guide_method_files), function(method) {
      x <- load_guide_results(cell_line, method)
      names(x)[names(x) == "p_value"] <- paste0(
        "p_", gsub("[^A-Za-z0-9]", "_", method)
      )
      x[, c(
        "gene", "guide",
        paste0("p_", gsub("[^A-Za-z0-9]", "_", method))
      )]
    })
    x <- Reduce(
      function(left, right) merge(left, right, by = c("gene", "guide")),
      method_tables
    )
    truth <- unique(guide_truth[
      guide_truth$cell_line == cell_line,
      c("gene", "truth")
    ])
    x <- merge(x, truth, by = "gene", all.x = TRUE)
    x$cell_line <- cell_line
    x
  }
))
competitor_p <- c("p_DESeq2", "p_edgeR_QL", "p_limma_voom")
proxy_null <- guide_disagreement$truth == 0L
forward <- proxy_null &
  guide_disagreement$p_BARCS >= 0.20 &
  apply(guide_disagreement[, competitor_p] < 0.05, 1L, all)
reverse <- proxy_null &
  guide_disagreement$p_BARCS < 0.05 &
  apply(guide_disagreement[, competitor_p] >= 0.05, 1L, all)
forward_candidates <- guide_disagreement[forward, , drop = FALSE]
forward_candidates$largest_competitor_p <- apply(
  forward_candidates[, competitor_p], 1L, max
)
best_within_cell_line <- do.call(rbind, lapply(
  split(forward_candidates, forward_candidates$cell_line),
  function(x) {
    x[order(x$largest_competitor_p, x$guide), , drop = FALSE][1L, ]
  }
))
best_within_cell_line$disagreement_ratio <-
  best_within_cell_line$p_BARCS /
  best_within_cell_line$largest_competitor_p
display_cases <- head(
  best_within_cell_line[
    order(best_within_cell_line$disagreement_ratio, decreasing = TRUE),
  ],
  3L
)
specified_key <- paste(
  case_specification$cell_line,
  case_specification$gene,
  case_specification$guide,
  sep = "::"
)
selected_key <- paste(
  display_cases$cell_line,
  display_cases$gene,
  display_cases$guide,
  sep = "::"
)
if (!setequal(specified_key, selected_key)) {
  stop("Displayed trajectories no longer match the documented selection rule.")
}
write.csv(
  data.frame(
    direction = c(
      "BARCS non-significant; all three alternatives significant",
      "BARCS significant; all three alternatives non-significant"
    ),
    selection_rule = c(
      "BARCS p >= 0.20; each alternative p < 0.05",
      "BARCS p < 0.05; each alternative p >= 0.05"
    ),
    proxy_null_guides = c(sum(forward, na.rm = TRUE), sum(reverse, na.rm = TRUE))
  ),
  file.path(
    "data", "derived", "liang_single_guide_disagreement_counts.csv"
  ),
  row.names = FALSE
)

load_count_trajectory <- function(cell_line, gene, guide) {
  cell_stem <- cell_stems[[cell_line]]
  path <- require_result_file(file.path(
    result_dir,
    paste0(cell_stem, "_longitudinal_counts.tsv")
  ))
  count_table <- read.delim(path, check.names = FALSE)
  count_columns <- setdiff(names(count_table), c("sgRNA", "Gene"))
  guide_row <- count_table[
    count_table$sgRNA == guide & count_table$Gene == gene,
    ,
    drop = FALSE
  ]
  if (nrow(guide_row) != 1L) {
    stop("Could not resolve one count row for ", cell_line, "/", guide)
  }

  library_sizes <- colSums(count_table[, count_columns, drop = FALSE])
  counts <- as.numeric(guide_row[, count_columns, drop = TRUE])
  samples <- count_columns
  data.frame(
    cell_line = cell_line,
    gene = gene,
    guide = guide,
    sample = samples,
    day = as.integer(sub(".*_Day([0-9]{2})_R[12]$", "\\1", samples)),
    replicate = as.integer(sub(".*_R([12])$", "\\1", samples)),
    processed_count = counts,
    library_size = library_sizes,
    counts_per_million = counts / library_sizes * 1e6,
    stringsAsFactors = FALSE
  )
}

count_trajectories <- do.call(rbind, lapply(
  seq_len(nrow(case_specification)),
  function(index) {
    load_count_trajectory(
      case_specification$cell_line[index],
      case_specification$gene[index],
      case_specification$guide[index]
    )
  }
))
write.csv(
  count_trajectories,
  file.path(
    "data", "derived", "liang_single_guide_count_trajectories.csv"
  ),
  row.names = FALSE
)

# Extend the original method-specific gene-level volcano plots across the four
# replicate-complete cell lines. Every method uses the same three-time-point
# longitudinal estimand with a replicate block within each cell line.
volcano_methods <- c(
  "BARCS", "MAGeCK-MLE", "edgeR-QL", "DESeq2", "limma-voom"
)
volcano_cell_lines <- names(cell_stems)
volcano_scores <- gene_scores[
  gene_scores$cell_line %in% volcano_cell_lines &
    gene_scores$method %in% volcano_methods &
    is.finite(gene_scores$effect) &
    is.finite(gene_scores$p_value) &
    is.finite(gene_scores$fdr),
  ,
  drop = FALSE
]
if (!identical(
  sort(unique(paste(
    volcano_scores$cell_line,
    volcano_scores$method,
    sep = "::"
  ))),
  sort(as.vector(outer(
    volcano_cell_lines,
    volcano_methods,
    paste,
    sep = "::"
  )))
)) {
  stop("All 20 cell-line-by-method longitudinal results are required.")
}

cell_line_colours <- c(
  HAP1 = "#0072B2",
  HEK293FT = "#009E73",
  `MDA-MB-231` = "#CC79A7",
  THP1 = "#D55E00"
)

pdf(
  file.path(
    "figures", "liang_longitudinal_volcano_trajectories.pdf"
  ),
  width = 10.5,
  height = 10.6,
  pointsize = 15,
  useDingbats = FALSE
)
layout(
  matrix(seq_len(9), nrow = 3, byrow = TRUE),
  heights = c(1, 1, 1.2)
)

# A signed logarithmic x transform expands the dense region around zero while
# retaining effect direction. The conventional -log10(p) axis is clipped so
# extreme values remain visible without defining the scale of every panel.
effect_scale <- 0.25
signed_log_effect <- function(value) {
  sign(value) * log10(1 + abs(value) / effect_scale)
}
y_clip <- 60
# Ticks are chosen so that their signed-log positions are close to evenly
# spaced. A denser set does not fit in a one-third-width panel, and base
# `axis()` then silently drops a different subset in each panel, which reads as
# five different x scales.
x_tick_effects <- c(-10, -2, -0.5, 0, 0.5, 2)
x_tick_positions <- signed_log_effect(x_tick_effects)
common_xlim <- range(signed_log_effect(volcano_scores$effect))
common_ylim <- c(0, y_clip)
# Name the transform on the axis rather than only in the caption, so the tick
# spacing is reproducible from the figure alone.
effect_axis_label <- sprintf(
  "Longitudinal effect (signed log, offset %s)",
  format(effect_scale)
)
# Every panel of Figure 1 uses colour for cell line and never for method, so
# the letters run straight through both rows.
volcano_letters <- LETTERS[seq_along(volcano_methods)]
trajectory_letters <- LETTERS[
  length(volcano_methods) + seq_len(nrow(case_specification))
]

# Sharing `common_xlim` is not sufficient for a shared x scale: `layout()` gives
# every cell the same width, so any panel that reserves extra width for a y-axis
# label ends up with a narrower plot region and therefore more effect units per
# inch. Keeping `mar` identical across all five panels is what makes the scale
# comparable; only the labels are conditional. `panel_plot_widths` asserts the
# invariant so the regression cannot return unnoticed.
volcano_mar <- c(4.5, 4.6, 2.8, 0.8)
panel_plot_widths <- numeric(length(volcano_methods))
panel_clipped_counts <- integer(length(volcano_methods))

for (method_index in seq_along(volcano_methods)) {
  method <- volcano_methods[method_index]
  panel <- volcano_scores[
    volcano_scores$method == method,
    ,
    drop = FALSE
  ]
  x_value <- signed_log_effect(panel$effect)
  raw_y_value <- -log10(pmax(panel$p_value, .Machine$double.xmin))
  y_value <- pmin(raw_y_value, y_clip)
  clipped <- raw_y_value > y_clip
  depleted_call <- panel$fdr < 0.10 & panel$effect < 0
  show_y_axis <- method_index %in% c(1L, 4L)

  par(
    mar = volcano_mar,
    cex.axis = 0.92,
    cex.lab = 0.98
  )
  plot(
    x_value,
    y_value,
    pch = 16,
    cex = 0.31,
    col = adjustcolor("#8C8C8C", alpha.f = 0.16),
    xlim = common_xlim,
    ylim = common_ylim,
    xaxt = "n",
    yaxt = if (show_y_axis) "s" else "n",
    xlab = effect_axis_label,
    ylab = if (show_y_axis) {
      expression(-log[10]("two-sided " * italic(p)))
    } else {
      ""
    },
    main = paste0(volcano_letters[method_index], "  ", method),
    bty = "l",
    cex.main = 1.05
  )
  panel_plot_widths[method_index] <- par("pin")[1L]
  panel_clipped_counts[method_index] <- sum(clipped)
  # A dense row of triangles on the top rule is a plotting artefact, not a pile
  # of identical p-values. Say how many points it stands for.
  mtext(
    sprintf("%d clipped", sum(clipped)),
    side = 3,
    line = 0.1,
    adj = 1,
    cex = 0.66,
    col = "#666666"
  )
  in_range <- x_tick_positions >= common_xlim[1L] &
    x_tick_positions <= common_xlim[2L]
  axis(
    1,
    at = x_tick_positions[in_range],
    labels = x_tick_effects[in_range]
  )
  if (!show_y_axis) {
    # Same breaks as the labelled panels, without repeating the numbers.
    axis(2, labels = FALSE, tcl = -0.3)
  }
  for (cell_line in volcano_cell_lines) {
    cell_rows <- depleted_call & panel$cell_line == cell_line
    points(
      x_value[cell_rows],
      y_value[cell_rows],
      pch = 16,
      cex = 0.38,
      col = adjustcolor(
        cell_line_colours[[cell_line]],
        alpha.f = 0.54
      )
    )
    clipped_rows <- cell_rows & clipped
    points(
      x_value[clipped_rows],
      y_value[clipped_rows],
      pch = 17,
      cex = 0.52,
      col = adjustcolor(
        cell_line_colours[[cell_line]],
        alpha.f = 0.70
      )
    )
  }
  points(
    x_value[clipped & !depleted_call],
    y_value[clipped & !depleted_call],
    pch = 17,
    cex = 0.48,
    col = adjustcolor("#777777", alpha.f = 0.55)
  )
  abline(h = -log10(0.05), col = "#3B5BDB", lty = 2)
  abline(v = 0, col = "#666666", lty = 3)
}

if (diff(range(panel_plot_widths)) > 1e-6) {
  stop(
    "Volcano panels have unequal plot widths, so their x axes are not on a ",
    "common scale: ",
    paste(sprintf("%s=%.3fin", volcano_methods, panel_plot_widths),
      collapse = ", "
    ),
    call. = FALSE
  )
}

write.csv(
  data.frame(
    method = volcano_methods,
    genes_plotted = as.integer(table(volcano_scores$method)[volcano_methods]),
    clipped_at_y = panel_clipped_counts,
    y_clip = y_clip,
    stringsAsFactors = FALSE
  ),
  file.path("data", "derived", "liang_volcano_clipped_counts.csv"),
  row.names = FALSE
)

par(mar = c(1.0, 1.0, 2.7, 1.0))
plot.new()
title(main = "Four-cell-line comparison", cex.main = 1.05)
legend(
  "center",
  title = "Colour: cell line, depleted at FDR 0.10",
  legend = c(
    "All genes",
    volcano_cell_lines,
    sprintf("Clipped at -log10(p) = %d", y_clip),
    "Nominal two-sided p = 0.05"
  ),
  pch = c(16, rep(16, length(volcano_cell_lines)), 17, NA),
  lty = c(NA, rep(NA, length(volcano_cell_lines)), NA, 2),
  col = c(
    "#8C8C8C",
    unname(cell_line_colours[volcano_cell_lines]),
    "#777777",
    "#3B5BDB"
  ),
  pt.cex = c(1.1, rep(1.1, length(volcano_cell_lines)), 1.0, NA),
  bty = "n",
  cex = 1.05
)

trajectory_shapes <- c(`1` = 21, `2` = 24)
trajectory_lines <- c(`1` = 1, `2` = 2)
trajectory_ylim <- c(0, 1.08 * max(count_trajectories$counts_per_million))
# Identical margins keep the three square panels the same size and aligned on a
# common left edge, for the same reason the volcano row needs it.
trajectory_mar <- c(4.6, 4.8, 3.9, 0.8)

for (case_index in seq_len(nrow(case_specification))) {
  cell_line <- case_specification$cell_line[case_index]
  guide <- case_specification$guide[case_index]
  trajectory <- count_trajectories[
    count_trajectories$cell_line == cell_line,
    ,
    drop = FALSE
  ]
  barcs_p <- case_results$p_value[
    case_results$cell_line == cell_line &
      case_results$method == "BARCS"
  ]
  # Filling by cell line keeps colour meaning the same thing in both rows. The
  # previous BARCS fill made one hue mean "cell line HAP1" above and "BARCS"
  # below, in a figure where the method is already given by the panel title.
  trajectory_fill <- cell_line_colours[[cell_line]]

  par(
    mar = trajectory_mar,
    cex.axis = 0.92,
    cex.lab = 0.98,
    pty = "s"
  )
  plot(
    NA,
    xlim = c(-0.8, 14.8),
    ylim = trajectory_ylim,
    xaxt = "n",
    xlab = "Day",
    ylab = if (case_index == 1L) {
      "Guide abundance (counts per million)"
    } else {
      ""
    },
    main = paste0(
      trajectory_letters[case_index], "  ",
      cell_line, ": ", guide
    ),
    bty = "l",
    cex.main = 1.05
  )
  axis(1, at = c(0, 7, 14))
  mtext(
    sprintf(
      "Two-sided p: BARCS %.2f (NS); others < 0.05",
      barcs_p
    ),
    side = 3,
    line = 0.15,
    cex = 0.76,
    col = "#444444"
  )

  for (replicate in sort(unique(trajectory$replicate))) {
    replicate_rows <- trajectory[
      trajectory$replicate == replicate,
      ,
      drop = FALSE
    ]
    replicate_rows <- replicate_rows[order(replicate_rows$day), ]
    lines(
      replicate_rows$day,
      replicate_rows$counts_per_million,
      lty = trajectory_lines[[as.character(replicate)]],
      lwd = 1.3,
      col = "#555555"
    )
    points(
      replicate_rows$day,
      replicate_rows$counts_per_million,
      pch = trajectory_shapes[[as.character(replicate)]],
      cex = 1.22,
      bg = trajectory_fill,
      col = "#222222",
      lwd = 0.8
    )
  }

  if (case_index == 1L) {
    legend(
      "topright",
      legend = c("Replicate 1", "Replicate 2"),
      pch = unname(trajectory_shapes),
      lty = unname(trajectory_lines),
      pt.bg = trajectory_fill,
      col = "#555555",
      pt.cex = 1.1,
      bty = "n",
      cex = 1.05
    )
  }
}
dev.off()

cat(
  "Figure 1 includes four-cell-line method volcano plots and three ",
  "guide-count trajectories.\n",
  sep = ""
)
