# Longitudinal HT-29 benchmark from Tzelepis et al., as used by Chronos.
#
# Raw counts are from Tzelepis Data S2 (PMC5081405/mmc6.zip). Published
# comparator effects, control lists, and the DepMap 20Q2 HT-29 expression
# profile are from Chronos Figshare article 14067047.  The analysis:
#   1. sums the three sequencing replicates at each late time point;
#   2. removes guides with fewer than 30 reads in pDNA, matching the published
#      Chronos preprocessing;
#   3. fits the same numeric time/25 design with beta-binomial regression and
#      the official MAGeCK-MLE executable;
#   4. compares both with the deposited Chronos-joint, MAGeCK-average, and
#      BAGEL2-average effects on the same gene/control universe; and
#   5. audits post-fit MAGeCK piecewise CNV correction using the HT-29
#      (ACH-000552) DepMap Public 20Q2 profile.

options(stringsAsFactors = FALSE)
source(file.path("R", "method_palette.R"))

benchmark_dir <- file.path("results", "chronos_tzelepis")
figure_path <- file.path("figures", "chronos_tzelepis_benchmark.pdf")
raw_dir <- file.path("data", "raw", "chronos")
dir.create(benchmark_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(figure_path), recursive = TRUE, showWarnings = FALSE)

required_files <- file.path(raw_dir, c(
  "Tzelepis_HT29_counts.txt.gz",
  "GeneFitnessEffect_ChronosJoint_Tzelepis.hdf5",
  "GeneFitnessEffect_MAGeCK_Tzelepis.hdf5",
  "GeneFitnessEffect_BAGEL2_Tzelepis.hdf5",
  "ReferenceEssentials.csv",
  "ReferenceNonEssentials.csv"
))
expression_path <- file.path(
  "data", "derived", "HT29_DepMap20Q2_expression.tsv"
)
required_files <- c(required_files, expression_path)
if (any(!file.exists(required_files))) {
  stop(
    "Chronos/Tzelepis source files are missing: ",
    paste(basename(required_files[!file.exists(required_files)]), collapse = ", "),
    call. = FALSE
  )
}
if (!requireNamespace("rhdf5", quietly = TRUE)) {
  stop("The Bioconductor package `rhdf5` is required.", call. = FALSE)
}

# Use the package's RcppArmadillo kernels when the local compiled library is
# present.  The base-R fallback remains valid but is slower.
cpp_library <- file.path("CB2", "src", paste0("CB2", .Platform$dynlib.ext))
if (file.exists(cpp_library) &&
    requireNamespace("Rcpp", quietly = TRUE) &&
    requireNamespace("RcppArmadillo", quietly = TRUE)) {
  suppressPackageStartupMessages(library(Rcpp))
  suppressPackageStartupMessages(library(RcppArmadillo))
  loaded_dlls <- names(getLoadedDLLs())
  if (!"CB2" %in% loaded_dlls) {
    dyn.load(cpp_library)
  }
  source(file.path("CB2", "R", "RcppExports.R"))
}
source(file.path("R", "bbreg.R"))

raw <- read.delim(
  gzfile(file.path(raw_dir, "Tzelepis_HT29_counts.txt.gz")),
  check.names = FALSE
)
days <- c(7, 10, 13, 16, 19, 22, 25)
late_counts <- sapply(days, function(day) {
  columns <- grep(
    sprintf("^HT29C3_d%d_rep", day), names(raw), value = TRUE
  )
  if (length(columns) != 3L) {
    stop("Expected three count columns for HT-29 day ", day, ".")
  }
  rowSums(raw[, columns, drop = FALSE])
})
counts <- cbind(pDNA = raw[["HumanV1-1"]], late_counts)
colnames(counts) <- c("pDNA", paste0("day", days))

# The three columns per day are sequencing/technical replicates of the same
# harvested population.  Summing prevents them from becoming false biological
# degrees of freedom in the Student t reference distribution.
keep <- counts[, "pDNA"] >= 30
counts <- counts[keep, , drop = FALSE]
guide <- raw$gRNA[keep]
gene <- raw$Gene[keep]
sample_data <- data.frame(
  sample = colnames(counts),
  day = c(0, days),
  time_25 = c(0, days / 25)
)

count_path <- file.path(benchmark_dir, "counts.tsv")
design_path <- file.path(benchmark_dir, "design.tsv")
write.table(
  data.frame(sgRNA = guide, Gene = gene, counts, check.names = FALSE),
  count_path, sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  data.frame(
    samples = sample_data$sample,
    baseline = 1,
    time_25 = sample_data$time_25
  ),
  design_path, sep = "\t", quote = FALSE, row.names = FALSE
)

bb_guide_path <- file.path(
  benchmark_dir, "beta_binomial_guide_results.csv.gz"
)
if (!file.exists(bb_guide_path) ||
    identical(Sys.getenv("RERUN_BB"), "1")) {
  start_time <- proc.time()
  bb_guide <- bb_screen(
    counts = counts,
    totals = colSums(counts),
    data = sample_data,
    formula = ~ time_25,
    term = "time_25",
    guide = guide,
    gene = gene,
    min_total_count = 30,
    ncores = 4
  )
  elapsed <- unname((proc.time() - start_time)[["elapsed"]])
  write.csv(bb_guide, gzfile(bb_guide_path), row.names = FALSE)
  write.csv(
    data.frame(
      n_guides = nrow(bb_guide),
      elapsed_seconds = elapsed,
      guides_per_second = nrow(bb_guide) / elapsed,
      converged = sum(bb_guide$converged, na.rm = TRUE)
    ),
    file.path(benchmark_dir, "beta_binomial_runtime.csv"),
    row.names = FALSE
  )
} else {
  bb_guide <- read.csv(gzfile(bb_guide_path))
}

combine_guides <- function(guide_result) {
  indices <- split(seq_len(nrow(guide_result)), guide_result$gene)
  result <- do.call(rbind, lapply(names(indices), function(gene_name) {
    index <- indices[[gene_name]]
    valid <- is.finite(guide_result$estimate[index]) &
      is.finite(guide_result$p_value[index])
    index <- index[valid]
    if (!length(index)) {
      return(NULL)
    }
    signed_z <- sign(guide_result$estimate[index]) *
      qnorm(
        pmax(
          guide_result$p_value[index] / 2,
          .Machine$double.xmin
        ),
        lower.tail = FALSE
      )
    combined_z <- sum(signed_z) / sqrt(length(signed_z))
    data.frame(
      gene = gene_name,
      n_guides = length(index),
      estimate = median(guide_result$estimate[index]),
      p_value = 2 * pnorm(-abs(combined_z)),
      converged_fraction = mean(guide_result$converged[index])
    )
  }))
  result$fdr <- p.adjust(result$p_value, method = "BH")
  rownames(result) <- NULL
  result
}

bb_gene <- combine_guides(bb_guide)
bb_gene_path <- file.path(
  benchmark_dir, "beta_binomial_gene_results.csv"
)
write.csv(bb_gene, bb_gene_path, row.names = FALSE)

source(file.path("R", "mageck.R"))
mageck <- mageck_executable()
mageck_py <- mageck_python()
mageck_check_version()
mageck_compat <- file.path("scripts", "mageck_compat.py")
mageck_cnv_correct <- file.path("scripts", "mageck_cnv_correct.py")
if (!all(file.exists(c(mageck_compat, mageck_cnv_correct)))) {
  stop("The MAGeCK helper scripts are missing from `scripts/`.", call. = FALSE)
}

cnv_path <- file.path("data", "derived", "HT29_DepMap20Q2_CNV.tsv")
if (!file.exists(cnv_path)) {
  stop("HT-29 copy-number profile is missing at `", cnv_path, "`.")
}

mageck_prefix <- file.path(benchmark_dir, "mageck_time")
mageck_cnv_prefix <- file.path(benchmark_dir, "mageck_time_cnv")
mageck_gene_path <- paste0(mageck_prefix, ".gene_summary.txt")
mageck_cnv_gene_path <- paste0(mageck_cnv_prefix, ".gene_summary.txt")
mageck_arguments <- c(
  "mle",
  "-k", count_path,
  "-d", design_path,
  "--norm-method", "median",
  "--permutation-round", "1",
  "--no-permutation-by-group",
  "--threads", "4"
)

if (!file.exists(mageck_gene_path) ||
    identical(Sys.getenv("RERUN_MAGECK"), "1")) {
  status <- system2(
    mageck,
    c(mageck_arguments, "-n", mageck_prefix)
  )
  if (status != 0 || !file.exists(mageck_gene_path)) {
    stop("Official continuous-time MAGeCK-MLE analysis failed.")
  }
}
if (!file.exists(mageck_cnv_gene_path) ||
    identical(Sys.getenv("RERUN_MAGECK"), "1")) {
  status <- system2(
    mageck_py,
    c(
      mageck_compat,
      mageck_arguments,
      "-n", mageck_cnv_prefix,
      "--cnv-norm", cnv_path,
      "--cell-line", "HT29_LARGE_INTESTINE"
    )
  )
  if (status != 0 || !file.exists(mageck_cnv_gene_path)) {
    stop("Official CNV-adjusted continuous-time MAGeCK-MLE analysis failed.")
  }
}

mageck_time <- read.delim(mageck_gene_path, check.names = FALSE)
mageck_time_cnv <- read.delim(mageck_cnv_gene_path, check.names = FALSE)

bb_effect_path <- file.path(
  benchmark_dir, "beta_binomial_gene_effects.tsv"
)
bb_cnv_effect_path <- file.path(
  benchmark_dir, "beta_binomial_gene_effects_cnv_corrected.tsv"
)
write.table(
  bb_gene[is.finite(bb_gene$estimate), c("gene", "estimate")],
  bb_effect_path, sep = "\t", quote = FALSE, row.names = FALSE
)
status <- system2(
  mageck_py,
  c(
    mageck_cnv_correct,
    "--effects", bb_effect_path,
    "--cnv", cnv_path,
    "--cell-line", "HT29_LARGE_INTESTINE",
    "--output", bb_cnv_effect_path
  )
)
if (status != 0 || !file.exists(bb_cnv_effect_path)) {
  stop("Official MAGeCK piecewise correction of beta-binomial effects failed.")
}
bb_cnv <- read.delim(bb_cnv_effect_path)

read_hdf5_effect <- function(filename, combine) {
  path <- file.path(raw_dir, filename)
  matrix <- rhdf5::h5read(path, "data")
  rows <- rhdf5::h5read(path, "dim_0")
  genes <- sub(
    " \\(.*$", "", rhdf5::h5read(path, "dim_1")
  )
  setNames(combine(matrix, rows), genes)
}

all_days_label <- "(7, 10, 13, 16, 19, 22, 25)"
chronos_joint <- read_hdf5_effect(
  "GeneFitnessEffect_ChronosJoint_Tzelepis.hdf5",
  function(matrix, rows) {
    matrix[, match(all_days_label, rows)]
  }
)
mageck_average <- read_hdf5_effect(
  "GeneFitnessEffect_MAGeCK_Tzelepis.hdf5",
  function(matrix, rows) {
    apply(matrix, 1, median, na.rm = TRUE)
  }
)
bagel_average <- read_hdf5_effect(
  "GeneFitnessEffect_BAGEL2_Tzelepis.hdf5",
  function(matrix, rows) {
    apply(matrix, 1, median, na.rm = TRUE)
  }
)

named_effect <- function(gene, effect) {
  setNames(as.numeric(effect), gene)
}
effects <- list(
  `BARCS time` = named_effect(bb_gene$gene, bb_gene$estimate),
  `Official MAGeCK time` = named_effect(
    mageck_time$Gene, mageck_time[["time_25|beta"]]
  ),
  `Chronos joint` = chronos_joint,
  `Published MAGeCK average` = mageck_average,
  `Published BAGEL2 average` = bagel_average
)
cnv_effects <- list(
  `BARCS time` = effects[["BARCS time"]],
  `BARCS time + CNV` = named_effect(
    bb_cnv$gene, bb_cnv$cnv_corrected_estimate
  ),
  `Official MAGeCK time` = effects[["Official MAGeCK time"]],
  `Official MAGeCK time + CNV` = named_effect(
    mageck_time_cnv$Gene,
    mageck_time_cnv[["time_25|beta"]]
  ),
  `Chronos joint` = chronos_joint,
  `Published MAGeCK average` = mageck_average,
  `Published BAGEL2 average` = bagel_average
)

read_reference <- function(filename) {
  values <- readLines(file.path(raw_dir, filename))
  # The deposited CSV has header "Gene" concatenated with its first value.
  values[1] <- sub("^Gene", "", values[1])
  sub(" \\(.*$", "", values)
}
essential <- read_reference("ReferenceEssentials.csv")
nonessential <- read_reference("ReferenceNonEssentials.csv")
expression <- read.delim(expression_path, check.names = FALSE)
unexpressed <- expression$SYMBOL[
  is.finite(expression$HT29_LARGE_INTESTINE) &
    expression$HT29_LARGE_INTESTINE < 0.5
]

# The Chronos time-course figure uses core essentials as positives and
# unexpressed HT-29 genes as its negative/null group.  Enforce one shared
# universe so no method benefits from missing difficult genes.
shared_genes <- Reduce(intersect, lapply(effects, names))
shared_genes <- intersect(
  shared_genes, union(essential, unexpressed)
)

curve_and_metrics <- function(effect) {
  effect <- effect[shared_genes]
  positive <- names(effect) %in% essential
  negative <- names(effect) %in% unexpressed & !positive
  valid <- is.finite(effect) & (positive | negative)
  effect <- effect[valid]
  positive <- positive[valid]
  negative <- negative[valid]

  order_index <- order(effect, na.last = NA)
  ordered_positive <- positive[order_index]
  true_positive <- cumsum(ordered_positive)
  false_positive <- cumsum(!ordered_positive)
  precision <- true_positive / seq_along(true_positive)
  recall <- true_positive / sum(positive)
  false_positive_rate <- false_positive / sum(negative)
  true_positive_rate <- recall

  ranks <- rank(effect, ties.method = "average")
  auc_positive_high <- (
    sum(ranks[positive]) -
      sum(positive) * (sum(positive) + 1) / 2
  ) / (sum(positive) * sum(negative))
  auroc <- 1 - auc_positive_high
  pr_auc <- sum(
    diff(c(0, recall)) *
      (head(c(1, precision), -1) + tail(c(1, precision), -1)) / 2
  )
  average_precision <- sum(
    true_positive[ordered_positive] /
      which(ordered_positive)
  ) / sum(positive)

  null <- effect[negative]
  true <- effect[positive]
  null_mad_about_mean <- mean(abs(null - mean(null)))
  nnmd <- (median(true) - mean(null)) / null_mad_about_mean
  cutoff <- unname(quantile(effect, 0.15, na.rm = TRUE))

  list(
    metrics = data.frame(
      n_essential = sum(positive),
      n_unexpressed = sum(negative),
      NNMD = nnmd,
      recall_at_90_precision = max(
        c(0, recall[precision >= 0.9])
      ),
      unexpressed_false_positives_15pct = sum(null < cutoff),
      AUROC = auroc,
      PR_AUC = pr_auc,
      average_precision = average_precision
    ),
    roc = data.frame(
      false_positive_rate = c(0, false_positive_rate, 1),
      true_positive_rate = c(0, true_positive_rate, 1)
    ),
    pr = data.frame(
      recall = c(0, recall),
      precision = c(1, precision)
    )
  )
}

evaluations <- lapply(effects, curve_and_metrics)
metrics <- do.call(rbind, lapply(names(evaluations), function(method) {
  cbind(method = method, evaluations[[method]]$metrics)
}))
rownames(metrics) <- NULL
write.csv(
  metrics, file.path(benchmark_dir, "benchmark_metrics.csv"),
  row.names = FALSE
)

effect_table <- data.frame(gene = shared_genes)
for (method in names(effects)) {
  effect_table[[method]] <- effects[[method]][shared_genes]
}
effect_table$is_reference_essential <- effect_table$gene %in% essential
effect_table$is_unexpressed_ht29 <- effect_table$gene %in% unexpressed
effect_table$is_reference_nonessential <-
  effect_table$gene %in% nonessential
write.csv(
  effect_table,
  gzfile(file.path(benchmark_dir, "shared_gene_effects.csv.gz")),
  row.names = FALSE
)

correlation <- cor(
  effect_table[, names(effects), drop = FALSE],
  method = "spearman", use = "pairwise.complete.obs"
)
correlation_table <- do.call(rbind, lapply(
  setdiff(names(effects), "BARCS time"),
  function(method) {
    data.frame(
      method_1 = "BARCS time",
      method_2 = method,
      n_genes = sum(
        is.finite(effect_table[["BARCS time"]]) &
          is.finite(effect_table[[method]])
      ),
      spearman = correlation["BARCS time", method]
    )
  }
))
write.csv(
  correlation_table,
  file.path(benchmark_dir, "effect_rank_correlations.csv"),
  row.names = FALSE
)

cnv <- read.delim(cnv_path)
cnv_bias <- do.call(rbind, lapply(names(cnv_effects), function(method) {
  effect <- cnv_effects[[method]]
  copy_number <- cnv$HT29_LARGE_INTESTINE[
    match(names(effect), cnv$SYMBOL)
  ]
  is_unexpressed <- names(effect) %in% unexpressed
  valid_all <- is.finite(effect) & is.finite(copy_number)
  valid_null <- valid_all & is_unexpressed
  data.frame(
    method = method,
    n_all = sum(valid_all),
    spearman_all = cor(
      effect[valid_all], copy_number[valid_all],
      method = "spearman"
    ),
    n_unexpressed = sum(valid_null),
    spearman_unexpressed = cor(
      effect[valid_null], copy_number[valid_null],
      method = "spearman"
    )
  )
}))
rownames(cnv_bias) <- NULL
write.csv(
  cnv_bias, file.path(benchmark_dir, "cnv_bias_diagnostic.csv"),
  row.names = FALSE
)

colors <- unname(barcs_method_colours[c(
  "BARCS", "MAGeCK", "Chronos", "Published MAGeCK", "BAGEL2"
)])
names(colors) <- names(effects)
line_types <- c(1, 1, 1, 2, 2)
names(line_types) <- names(effects)

pdf(figure_path, width = 10, height = 8)
old_par <- par(
  mfrow = c(2, 2), mar = c(4.2, 4.4, 3.0, 1.0),
  las = 1, xaxs = "i", yaxs = "i"
)

plot(
  0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
  xlab = "False-positive rate", ylab = "True-positive rate",
  main = "(A) Reference-essential discrimination"
)
abline(0, 1, col = "grey80", lty = 3)
for (method in names(evaluations)) {
  lines(
    evaluations[[method]]$roc,
    col = colors[[method]], lwd = 2, lty = line_types[[method]]
  )
}
legend(
  "bottomright", legend = names(effects), col = colors,
  lty = line_types, lwd = 2, cex = 0.72, bty = "n"
)

plot(
  0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
  xlab = "Recall", ylab = "Precision",
  main = "(B) Essential vs unexpressed genes"
)
abline(h = metrics$n_essential[1] /
  (metrics$n_essential[1] + metrics$n_unexpressed[1]),
col = "grey80", lty = 3)
for (method in names(evaluations)) {
  lines(
    evaluations[[method]]$pr,
    col = colors[[method]], lwd = 2, lty = line_types[[method]]
  )
}

barplot(
  metrics$recall_at_90_precision,
  names.arg = c("BARCS", "MAGeCK\ntime", "Chronos\njoint",
                "MAGeCK\naverage", "BAGEL2\naverage"),
  col = unname(colors[metrics$method]), border = NA,
  ylim = c(0, 1), ylab = "Recall at >= 90% precision",
  main = "(C) High-precision recovery", las = 2, cex.names = 0.72
)
abline(h = seq(0, 1, 0.2), col = "white", lwd = 0.7)

cnv_plot_methods <- c(
  "BARCS time", "BARCS time + CNV",
  "Official MAGeCK time", "Official MAGeCK time + CNV",
  "Chronos joint", "Published MAGeCK average",
  "Published BAGEL2 average"
)
cnv_values <- abs(cnv_bias$spearman_unexpressed[
  match(cnv_plot_methods, cnv_bias$method)
])
cnv_labels <- c(
  "BARCS", "BARCS + CNV", "MAGeCK time", "MAGeCK time + CNV",
  "Chronos joint", "MAGeCK average", "BAGEL2 average"
)
cnv_colors <- c(
  colors[["BARCS time"]], colors[["BARCS time"]],
  colors[["Official MAGeCK time"]], colors[["Official MAGeCK time"]],
  colors[["Chronos joint"]], colors[["Published MAGeCK average"]],
  colors[["Published BAGEL2 average"]]
)
cnv_pch <- c(16, 1, 16, 1, 16, 16, 16)
plot(
  0, 0, type = "n",
  xlim = c(0, max(cnv_values) * 1.08),
  ylim = c(0.5, length(cnv_values) + 0.5),
  xlab = "|Spearman(effect, copy number)| (closer to 0 is better)",
  ylab = "", yaxt = "n",
  main = "(D) Residual CNV association"
)
axis(2, at = seq_along(cnv_values), labels = rev(cnv_labels),
     las = 1, cex.axis = 0.68)
abline(v = pretty(c(0, cnv_values)), col = "grey92", lwd = 0.8)
segments(
  x0 = 0, y0 = seq_along(cnv_values),
  x1 = rev(cnv_values), y1 = seq_along(cnv_values),
  col = rev(cnv_colors), lwd = 2
)
points(
  rev(cnv_values), seq_along(cnv_values),
  col = rev(cnv_colors), pch = rev(cnv_pch), cex = 1.15, lwd = 2
)

par(old_par)
dev.off()

cat("\nChronos/Tzelepis longitudinal benchmark\n")
cat(
  nrow(counts), "guides,", length(unique(gene)), "genes,",
  ncol(counts), "aggregated time points\n"
)
print(metrics, row.names = FALSE)
cat("\nEffect-rank correlation with BARCS time slope\n")
print(correlation_table, row.names = FALSE)
cat("\nCNV diagnostic\n")
print(cnv_bias, row.names = FALSE)
