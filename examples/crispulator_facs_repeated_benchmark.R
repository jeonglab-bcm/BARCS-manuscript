#!/usr/bin/env Rscript

# Repeat the Crispulator FACS benchmark across prespecified simulation seeds.
#
# The per-seed count matrices and MAGeCK intermediates stay under `results/`
# (git-ignored). The compact cross-seed metrics and manuscript figure are
# written to versioned paths.

options(stringsAsFactors = FALSE)
source(file.path("R", "method_palette.R"))

n_seeds <- as.integer(Sys.getenv("CRISPULATOR_N_SEEDS", "5"))
if (!is.finite(n_seeds) || n_seeds < 2L) {
  stop("`CRISPULATOR_N_SEEDS` must be an integer of at least two.")
}
seeds <- 20250724L + seq.int(0L, n_seeds - 1L)
moi <- as.numeric(Sys.getenv("CRISPULATOR_MOI", "0.25"))
high_quality_fraction <- as.numeric(Sys.getenv(
  "CRISPULATOR_HIGH_QUALITY_GUIDE_FRACTION", "0.90"
))
if (!is.finite(moi) || moi <= 0 || moi >= 0.5) {
  stop("`CRISPULATOR_MOI` must be greater than zero and below 0.5.")
}
if (!is.finite(high_quality_fraction) ||
    high_quality_fraction < 0 || high_quality_fraction > 1) {
  stop(
    "`CRISPULATOR_HIGH_QUALITY_GUIDE_FRACTION` must be between zero and one."
  )
}
scenario <- sprintf(
  "moi_%s_quality_%s",
  format(moi, trim = TRUE, scientific = FALSE),
  format(high_quality_fraction, trim = TRUE, scientific = FALSE)
)
root <- file.path("results", "crispulator_facs", "repeated", scenario)
dir.create(root, recursive = TRUE, showWarnings = FALSE)

all_metrics <- lapply(seeds, function(seed) {
  run_dir <- file.path(root, paste0("seed_", seed))
  data_dir <- file.path(run_dir, "data")
  result_dir <- file.path(run_dir, "analysis")
  figure_path <- file.path(run_dir, "benchmark.pdf")
  metric_file <- file.path(result_dir, "benchmark_metrics.csv")
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(metric_file) ||
      identical(Sys.getenv("RERUN_CRISPULATOR"), "1")) {
    simulation_status <- system2(
      "julia",
      c(
        "--project=julia",
        file.path("julia", "simulate_crispulator_facs.jl"),
        data_dir,
        "4",
        as.character(seed),
        as.character(moi),
        as.character(high_quality_fraction)
      ),
      stdout = FALSE,
      stderr = FALSE
    )
    if (simulation_status != 0) {
      stop("Crispulator simulation failed for seed ", seed, ".")
    }

    analysis_status <- system2(
      "Rscript",
      file.path("examples", "crispulator_facs_benchmark.R"),
      env = c(
        paste0("CRISPULATOR_DATA_DIR=", data_dir),
        paste0("CRISPULATOR_RESULT_DIR=", result_dir),
        paste0("CRISPULATOR_FIGURE_PATH=", figure_path)
      ),
      stdout = FALSE,
      stderr = FALSE
    )
    if (analysis_status != 0) {
      stop("BARCS/MAGeCK analysis failed for seed ", seed, ".")
    }
  }
  metrics <- read.csv(metric_file)
  metrics$seed <- seed
  metrics$moi <- moi
  metrics$high_quality_guide_fraction <- high_quality_fraction
  metrics
})
all_metrics <- do.call(rbind, all_metrics)

design_concordance <- do.call(rbind, lapply(seeds, function(seed) {
  result_dir <- file.path(
    root, paste0("seed_", seed), "analysis"
  )
  comparisons <- list(
    BARCS = c(
      "barcs_low_bulk_high_gene_results.csv",
      "barcs_two_tails_gene_results.csv"
    ),
    `MAGeCK-MLE` = c(
      "mageck_mle_low_bulk_high_gene_results.csv",
      "mageck_mle_two_tails_gene_results.csv"
    )
  )
  do.call(rbind, lapply(names(comparisons), function(method) {
    paths <- file.path(result_dir, comparisons[[method]])
    three_sample <- read.csv(paths[1L])
    tails <- read.csv(paths[2L])
    shared <- merge(
      three_sample[, c("gene", "estimate", "p_value")],
      tails[, c("gene", "estimate", "p_value")],
      by = "gene",
      suffixes = c("_three", "_tails")
    )
    data.frame(
      seed = seed,
      method = method,
      moi = moi,
      high_quality_guide_fraction = high_quality_fraction,
      effect_spearman = cor(
        shared$estimate_three,
        shared$estimate_tails,
        method = "spearman"
      ),
      significance_spearman = cor(
        -log10(pmax(shared$p_value_three, .Machine$double.xmin)),
        -log10(pmax(shared$p_value_tails, .Machine$double.xmin)),
        method = "spearman"
      )
    )
  }))
}))
write.csv(
  design_concordance,
  file.path(
    "data", "derived", "crispulator_facs_design_concordance.csv"
  ),
  row.names = FALSE
)

metric_path <- file.path(
  "data", "derived", "crispulator_facs_repeated_metrics.csv"
)
write.csv(all_metrics, metric_path, row.names = FALSE)

metric_columns <- c(
  "auroc",
  "average_precision",
  "effect_spearman_active",
  "direction_accuracy_active",
  "directional_recall_fdr_0_10",
  "empirical_fdp_fdr_0_10",
  "f1_fdr_0_10",
  "negative_control_p_below_0_05"
)
summary_rows <- do.call(rbind, lapply(
  split(all_metrics, interaction(
    all_metrics$method, all_metrics$design, drop = TRUE
  )),
  function(group) {
    result <- data.frame(
      method = group$method[1L],
      design = group$design[1L],
      simulations = nrow(group),
      moi = moi,
      high_quality_guide_fraction = high_quality_fraction
    )
    for (metric in metric_columns) {
      result[[paste0(metric, "_mean")]] <- mean(group[[metric]])
      result[[paste0(metric, "_sd")]] <- sd(group[[metric]])
    }
    result
  }
))
rownames(summary_rows) <- NULL
write.csv(
  summary_rows,
  file.path(
    "data", "derived", "crispulator_facs_repeated_summary.csv"
  ),
  row.names = FALSE
)

figure_path <- file.path("figures", "crispulator_facs_benchmark.pdf")
pdf(figure_path, width = 9, height = 3.6, useDingbats = FALSE)
layout(matrix(1:3, nrow = 1), widths = c(1.08, 1.25, 1))

# A: requested low/bulk/high sampling geometry.
par(mar = c(4.2, 4.2, 2.3, 0.8))
q25 <- qnorm(0.25)
x <- seq(-3.5, 3.5, length.out = 800)
y <- dnorm(x)
plot(
  x, y, type = "n", xlab = "Latent FACS phenotype (z)",
  ylab = "Cell density", yaxt = "n", bty = "l",
  main = "A  Sampling design"
)
regions <- list(c(-3.5, q25), c(-q25, 3.5))
fills <- c("#56B4E980", "#D55E0080")
for (index in seq_along(regions)) {
  keep <- x >= regions[[index]][1L] & x <= regions[[index]][2L]
  polygon(
    c(x[keep], rev(x[keep])),
    c(y[keep], rep(0, sum(keep))),
    col = fills[index], border = NA
  )
}
lines(x, y, lwd = 1.5)
abline(v = c(q25, -q25), lty = 3, col = "#555555")
text(
  c(-1.45, 1.45), c(0.07, 0.07),
  c("Low\n25%", "High\n25%"),
  cex = 0.9
)
arrows(-3.25, 0.31, 3.25, 0.31, code = 3, length = 0.05)
text(0, 0.335, "Bulk 0-100% (overlapping)", cex = 0.84)

# B: cross-seed average precision and directional recall.
sorted <- all_metrics[all_metrics$method != "Bulk vs input", ]
sorted$label <- paste(sorted$method, sorted$design, sep = "\n")
label_order <- c(
  "BARCS\nLow + bulk + high",
  "MAGeCK-MLE\nLow + bulk + high",
  "BARCS\nTwo 25% tails",
  "MAGeCK-MLE\nTwo 25% tails"
)
display_labels <- c(
  "BARCS\n3 samples",
  "MAGeCK-MLE\n3 samples",
  "BARCS\n2 tails",
  "MAGeCK-MLE\n2 tails"
)
sorted$label <- factor(sorted$label, levels = label_order)
means <- do.call(rbind, lapply(split(sorted, sorted$label), function(group) {
  data.frame(
    label = group$label[1L],
    ap_mean = mean(group$average_precision),
    ap_sd = sd(group$average_precision),
    recall_mean = mean(group$directional_recall_fdr_0_10),
    recall_sd = sd(group$directional_recall_fdr_0_10)
  )
}))
means <- means[match(label_order, as.character(means$label)), ]
values <- rbind(
  `Average precision` = means$ap_mean,
  `Directional recall` = means$recall_mean
)
par(mar = c(7.5, 4.2, 2.3, 0.8))
positions <- barplot(
  values,
  beside = TRUE,
  col = rep(
    c(
      barcs_method_colours[["BARCS"]],
      barcs_method_colours[["MAGeCK"]],
      barcs_method_colours[["BARCS"]],
      barcs_method_colours[["MAGeCK"]]
    ),
    each = 2L
  ),
  density = rep(c(-1, 35), times = 4L),
  angle = 45,
  border = rep(
    c(
      barcs_method_colours[["BARCS"]],
      barcs_method_colours[["MAGeCK"]],
      barcs_method_colours[["BARCS"]],
      barcs_method_colours[["MAGeCK"]]
    ),
    each = 2L
  ),
  ylim = c(0, 1.15),
  ylab = "Truth-recovery metric",
  names.arg = display_labels,
  las = 2,
  cex.names = 0.72,
  main = sprintf("B  Recovery across %d seeds", n_seeds)
)
error_sd <- rbind(means$ap_sd, means$recall_sd)
has_error <- is.finite(error_sd) & error_sd > 0
arrows(
  positions[has_error],
  (values - error_sd)[has_error],
  positions[has_error],
  (values + error_sd)[has_error],
  angle = 90,
  code = 3,
  length = 0.035
)
legend(
  "top",
  legend = rownames(values),
  fill = c("#555555", "#555555"),
  density = c(-1, 35),
  angle = 45,
  border = "#555555",
  horiz = TRUE,
  bty = "n",
  cex = 0.7
)

# C: the 0--100% bulk sample carries no directional FACS signal.
bulk <- all_metrics[all_metrics$method == "Bulk vs input", ]
par(mar = c(5.2, 4.2, 2.3, 0.8))
bulk_values <- cbind(
  AUROC = bulk$auroc,
  Spearman = bulk$effect_spearman_active
)
boxplot(
  bulk_values,
  col = c("#99999966", "#99999966"),
  border = "#555555",
  ylim = c(-0.2, 0.7),
  ylab = "Null-reference metric",
  main = "C  Bulk versus input",
  names = c("Active-gene\nAUROC", "Effect\nSpearman")
)
abline(h = c(0, 0.5), lty = 3, col = "#777777")
stripchart(
  bulk_values,
  vertical = TRUE,
  method = "jitter",
  add = TRUE,
  pch = 16,
  cex = 0.65,
  col = "#333333"
)
text(
  c(1, 2),
  c(0.65, 0.65),
  sprintf(
    "mean %.2f",
    c(mean(bulk$auroc), mean(bulk$effect_spearman_active))
  ),
  cex = 0.72
)
dev.off()

print(summary_rows)
