#!/usr/bin/env Rscript

# Manuscript Figure 1: longitudinal Liang method comparison.
#
# The upper rows show the four replicate-complete cell lines in method-specific
# volcano plots for
# BARCS, MAGeCK-MLE, edgeR-QL, DESeq2, and limma-voom. Panels A-C show the
# observed day 0, 7, and 14 guide-abundance trajectories for representative
# false-positive calls.

options(stringsAsFactors = FALSE)
source(file.path("R", "method_palette.R"))
source(file.path("R", "bbreg.R"))

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
  cell_line = c("HAP1", "THP1", "MDA-MB-231"),
  gene = c(
    "Hum_XLOC_047119",
    "Hum_XLOC_037359",
    "Hum_XLOC_047245"
  ),
  guide = c("gL_028460", "gL_044332", "gL_028545"),
  stringsAsFactors = FALSE
)

gene_scores <- read.csv(gzfile(score_path))
guide_truth <- unique(gene_scores[
  gene_scores$method == "BARCS",
  c("cell_line", "gene", "truth")
])

load_guide_results <- function(cell_line, method) {
  cell_stem <- cell_stems[[cell_line]]
  path <- file.path(
    result_dir,
    paste0(
      cell_stem, "_longitudinal_",
      guide_method_files[[method]], "_guide.csv.gz"
    )
  )
  result <- read.csv(gzfile(path))
  if (method == "BARCS") {
    controls <- readLines(file.path(
      result_dir,
      paste0(cell_stem, "_longitudinal_non_targeting_guides.txt")
    ))
    result <- bb_calibrate_controls(result, result$guide %in% controls)
  }
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

load_count_trajectory <- function(cell_line, gene, guide) {
  cell_stem <- cell_stems[[cell_line]]
  path <- file.path(
    result_dir,
    paste0(cell_stem, "_longitudinal_counts.tsv")
  )
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

barcs_colour <- barcs_method_colours[["BARCS"]]
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
  height = 9.4,
  pointsize = 15,
  useDingbats = FALSE
)
layout(
  matrix(seq_len(9), nrow = 3, byrow = TRUE),
  heights = c(1, 1, 0.95)
)

# A signed logarithmic x transform expands the dense region around zero while
# retaining effect direction. A second logarithm on the usual -log10(p) axis
# keeps extreme values visible without allowing them to dominate the panels.
effect_scale <- 0.25
signed_log_effect <- function(value) {
  sign(value) * log10(1 + abs(value) / effect_scale)
}
double_log_p <- function(value) {
  log10(1 - log10(pmax(value, .Machine$double.xmin)))
}
x_tick_effects <- c(-10, -5, -2, -1, -0.5, 0, 0.5, 1, 2)
x_tick_positions <- signed_log_effect(x_tick_effects)
y_tick_logp <- c(0, 1, 2, 5, 10, 20, 60, 300)
y_tick_positions <- log10(1 + y_tick_logp)
common_xlim <- range(signed_log_effect(volcano_scores$effect))
common_ylim <- range(double_log_p(volcano_scores$p_value))

for (method_index in seq_along(volcano_methods)) {
  method <- volcano_methods[method_index]
  panel <- volcano_scores[
    volcano_scores$method == method,
    ,
    drop = FALSE
  ]
  x_value <- signed_log_effect(panel$effect)
  y_value <- double_log_p(panel$p_value)
  depleted_call <- panel$fdr < 0.10 & panel$effect < 0
  show_y_axis <- method_index %in% c(1L, 4L)

  par(
    mar = c(4.5, if (show_y_axis) 4.6 else 1.2, 2.8, 0.8),
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
    yaxt = "n",
    xlab = "Longitudinal time-slope effect (signed log scale)",
    ylab = if (show_y_axis) {
      expression(log[10](1-log[10]("two-sided " * italic(p))))
    } else {
      ""
    },
    main = method,
    bty = "l",
    cex.main = 1.05
  )
  axis(
    1,
    at = x_tick_positions[
      x_tick_positions >= common_xlim[1L] &
        x_tick_positions <= common_xlim[2L]
    ],
    labels = x_tick_effects[
      x_tick_positions >= common_xlim[1L] &
        x_tick_positions <= common_xlim[2L]
    ]
  )
  if (show_y_axis) {
    axis(
      2,
      at = y_tick_positions[
        y_tick_positions >= common_ylim[1L] &
          y_tick_positions <= common_ylim[2L]
      ],
      labels = y_tick_logp[
        y_tick_positions >= common_ylim[1L] &
          y_tick_positions <= common_ylim[2L]
      ]
    )
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
  }
  abline(h = double_log_p(0.05), col = "#3B5BDB", lty = 2)
  abline(v = 0, col = "#666666", lty = 3)
}

par(mar = c(1.0, 1.0, 2.7, 1.0))
plot.new()
title(main = "Four-cell-line comparison", cex.main = 1.05)
legend(
  "center",
  legend = c(
    "All genes",
    volcano_cell_lines,
    "Nominal two-sided p = 0.05"
  ),
  pch = c(16, rep(16, length(volcano_cell_lines)), NA),
  lty = c(NA, rep(NA, length(volcano_cell_lines)), 2),
  col = c(
    "#8C8C8C",
    unname(cell_line_colours[volcano_cell_lines]),
    "#3B5BDB"
  ),
  pt.cex = c(0.8, rep(0.8, length(volcano_cell_lines)), NA),
  bty = "n",
  cex = 0.78
)

trajectory_shapes <- c(`1` = 21, `2` = 24)
trajectory_lines <- c(`1` = 1, `2` = 2)
trajectory_ylim <- c(0, 1.08 * max(count_trajectories$counts_per_million))

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

  par(
    mar = c(4.6, if (case_index == 1L) 4.8 else 2.2, 3.9, 0.8),
    cex.axis = 0.94,
    cex.lab = 1.0
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
      LETTERS[case_index], "  ",
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
      bg = barcs_colour,
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
      pt.bg = barcs_colour,
      col = "#555555",
      pt.cex = 0.9,
      bty = "n",
      cex = 0.82
    )
  }
}
dev.off()

cat(
  "Figure 1 includes four-cell-line method volcano plots and three ",
  "guide-count trajectories.\n",
  sep = ""
)
