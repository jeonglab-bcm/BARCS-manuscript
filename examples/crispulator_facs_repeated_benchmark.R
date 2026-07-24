#!/usr/bin/env Rscript

# Repeat the Crispulator FACS benchmark across a parameter grid and
# prespecified simulation seeds.
#
# The per-seed count matrices and MAGeCK intermediates stay under `results/`
# (git-ignored). Compact cross-seed metrics and figures are versioned. The
# default one-at-a-time grid varies MOI, high-quality-guide fraction, and gene
# count around the manuscript baseline without confounding their effects.

options(stringsAsFactors = FALSE)
source(file.path("R", "method_palette.R"))

parse_numeric_values <- function(name, default) {
  text <- Sys.getenv(name, paste(default, collapse = ","))
  values <- as.numeric(trimws(strsplit(text, ",", fixed = TRUE)[[1L]]))
  if (!length(values) || any(!is.finite(values))) {
    stop("`", name, "` must be a comma-separated numeric list.")
  }
  unique(values)
}

format_parameter <- function(value) {
  format(value, trim = TRUE, scientific = FALSE)
}

cache_scenario_name <- function(moi, quality, genes) {
  name <- sprintf(
    "moi_%s_quality_%s",
    format_parameter(moi),
    format_parameter(quality)
  )
  if (genes != 400L) {
    name <- paste0(name, "_genes_", genes)
  }
  name
}

n_seeds <- as.integer(Sys.getenv("CRISPULATOR_N_SEEDS", "5"))
if (!is.finite(n_seeds) || n_seeds < 2L) {
  stop("`CRISPULATOR_N_SEEDS` must be an integer of at least two.")
}
seeds <- 20250724L + seq.int(0L, n_seeds - 1L)
baseline_moi <- as.numeric(Sys.getenv("CRISPULATOR_MOI", "0.25"))
baseline_quality <- as.numeric(Sys.getenv(
  "CRISPULATOR_HIGH_QUALITY_GUIDE_FRACTION", "0.90"
))
baseline_genes <- as.integer(Sys.getenv("CRISPULATOR_GENES", "400"))
moi_values <- parse_numeric_values(
  "CRISPULATOR_MOI_VALUES", c(0.10, 0.25, 0.40)
)
quality_values <- parse_numeric_values(
  "CRISPULATOR_HIGH_QUALITY_GUIDE_FRACTION_VALUES", c(0.60, 0.75, 0.90)
)
gene_values <- as.integer(parse_numeric_values(
  "CRISPULATOR_GENE_VALUES", c(100, 400, 1000)
))
grid_mode <- Sys.getenv("CRISPULATOR_GRID_MODE", "one_at_a_time")

if (!is.finite(baseline_moi) ||
    baseline_moi <= 0 || baseline_moi >= 0.5 ||
    any(moi_values <= 0 | moi_values >= 0.5)) {
  stop("`CRISPULATOR_MOI` must be greater than zero and below 0.5.")
}
if (!is.finite(baseline_quality) ||
    baseline_quality < 0 || baseline_quality > 1 ||
    any(quality_values < 0 | quality_values > 1)) {
  stop(
    "Guide-quality fractions must be between zero and one."
  )
}
if (!is.finite(baseline_genes) || baseline_genes < 20L ||
    any(gene_values < 20L)) {
  stop("Gene counts must be integers of at least 20.")
}
if (!grid_mode %in% c("single", "one_at_a_time", "full_factorial")) {
  stop(
    "`CRISPULATOR_GRID_MODE` must be single, one_at_a_time, or full_factorial."
  )
}

baseline <- data.frame(
  moi = baseline_moi,
  high_quality_guide_fraction = baseline_quality,
  genes = baseline_genes
)
scenarios <- switch(
  grid_mode,
  single = baseline,
  one_at_a_time = unique(rbind(
    baseline,
    data.frame(
      moi = moi_values,
      high_quality_guide_fraction = baseline_quality,
      genes = baseline_genes
    ),
    data.frame(
      moi = baseline_moi,
      high_quality_guide_fraction = quality_values,
      genes = baseline_genes
    ),
    data.frame(
      moi = baseline_moi,
      high_quality_guide_fraction = baseline_quality,
      genes = gene_values
    )
  )),
  full_factorial = expand.grid(
    moi = moi_values,
    high_quality_guide_fraction = quality_values,
    genes = gene_values,
    KEEP.OUT.ATTRS = FALSE
  )
)
scenarios$scenario_id <- sprintf(
  "moi_%s_quality_%s_genes_%d",
  vapply(scenarios$moi, format_parameter, character(1L)),
  vapply(
    scenarios$high_quality_guide_fraction,
    format_parameter,
    character(1L)
  ),
  scenarios$genes
)
scenarios$is_baseline <- with(
  scenarios,
  moi == baseline_moi &
    high_quality_guide_fraction == baseline_quality &
    genes == baseline_genes
)
scenarios <- scenarios[
  order(
    !scenarios$is_baseline,
    scenarios$genes,
    scenarios$high_quality_guide_fraction,
    scenarios$moi
  ),
]
write.csv(
  scenarios,
  file.path("data", "derived", "crispulator_facs_parameter_grid.csv"),
  row.names = FALSE
)

root <- file.path("results", "crispulator_facs", "repeated")
dir.create(root, recursive = TRUE, showWarnings = FALSE)

runs <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(index) {
  cbind(
    scenarios[rep(index, length(seeds)), ],
    seed = seeds,
    row.names = NULL
  )
}))

cache_matches <- function(parameter_file, run) {
  if (!file.exists(parameter_file)) {
    return(FALSE)
  }
  parameters <- read.delim(parameter_file)
  values <- setNames(parameters$value, parameters$parameter)
  isTRUE(all.equal(
    as.numeric(values[["multiplicity_of_infection"]]),
    run$moi
  )) &&
    isTRUE(all.equal(
      as.numeric(values[["high_quality_guide_fraction"]]),
      run$high_quality_guide_fraction
    )) &&
    identical(as.integer(values[["genes"]]), as.integer(run$genes)) &&
    identical(as.integer(values[["seed"]]), as.integer(run$seed))
}

all_metrics <- lapply(seq_len(nrow(runs)), function(index) {
  run <- runs[index, ]
  # Reuse the earlier cache layout for 400-gene scenarios after validating
  # parameters.tsv; other gene counts receive an explicit suffix.
  cache_scenario <- cache_scenario_name(
    run$moi,
    run$high_quality_guide_fraction,
    run$genes
  )
  run_dir <- file.path(
    root, cache_scenario, paste0("seed_", run$seed)
  )
  data_dir <- file.path(run_dir, "data")
  result_dir <- file.path(run_dir, "analysis")
  figure_path <- file.path(run_dir, "benchmark.pdf")
  metric_file <- file.path(result_dir, "benchmark_metrics.csv")
  parameter_file <- file.path(data_dir, "parameters.tsv")
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(metric_file) ||
      !cache_matches(parameter_file, run) ||
      identical(Sys.getenv("RERUN_CRISPULATOR"), "1")) {
    message(
      sprintf(
        "[%d/%d] MOI=%s, quality=%s, genes=%d, seed=%d",
        index,
        nrow(runs),
        format_parameter(run$moi),
        format_parameter(run$high_quality_guide_fraction),
        run$genes,
        run$seed
      )
    )
    simulation_status <- system2(
      "julia",
      c(
        "--project=julia",
        file.path("julia", "simulate_crispulator_facs.jl"),
        data_dir,
        "4",
        as.character(run$seed),
        as.character(run$moi),
        as.character(run$high_quality_guide_fraction),
        as.character(run$genes)
      ),
      stdout = FALSE,
      stderr = FALSE
    )
    if (simulation_status != 0) {
      stop("Crispulator simulation failed for ", run$scenario_id, ".")
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
      stop("BARCS/MAGeCK analysis failed for ", run$scenario_id, ".")
    }
  }
  metrics <- read.csv(metric_file)
  metrics$scenario_id <- run$scenario_id
  metrics$seed <- run$seed
  metrics$moi <- run$moi
  metrics$high_quality_guide_fraction <- run$high_quality_guide_fraction
  metrics$genes <- run$genes
  metrics$is_baseline <- run$is_baseline
  metrics
})
all_metrics <- do.call(rbind, all_metrics)

design_concordance <- do.call(rbind, lapply(seq_len(nrow(runs)), function(
    index) {
  run <- runs[index, ]
  cache_scenario <- cache_scenario_name(
    run$moi,
    run$high_quality_guide_fraction,
    run$genes
  )
  result_dir <- file.path(
    root, cache_scenario, paste0("seed_", run$seed), "analysis"
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
      scenario_id = run$scenario_id,
      seed = run$seed,
      method = method,
      moi = run$moi,
      high_quality_guide_fraction = run$high_quality_guide_fraction,
      genes = run$genes,
      is_baseline = run$is_baseline,
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
    "data", "derived", "crispulator_facs_parameter_grid_concordance.csv"
  ),
  row.names = FALSE
)
write.csv(
  all_metrics,
  file.path(
    "data", "derived", "crispulator_facs_parameter_grid_metrics.csv"
  ),
  row.names = FALSE
)

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
    all_metrics$scenario_id,
    all_metrics$method,
    all_metrics$design,
    drop = TRUE
  )),
  function(group) {
    result <- data.frame(
      scenario_id = group$scenario_id[1L],
      method = group$method[1L],
      design = group$design[1L],
      simulations = nrow(group),
      moi = group$moi[1L],
      high_quality_guide_fraction =
        group$high_quality_guide_fraction[1L],
      genes = group$genes[1L],
      is_baseline = group$is_baseline[1L]
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
    "data", "derived", "crispulator_facs_parameter_grid_summary.csv"
  ),
  row.names = FALSE
)

# Preserve the established manuscript files as the baseline-only subset.
baseline_metrics <- all_metrics[all_metrics$is_baseline, ]
baseline_concordance <- design_concordance[design_concordance$is_baseline, ]
baseline_summary <- summary_rows[summary_rows$is_baseline, ]
write.csv(
  baseline_metrics,
  file.path(
    "data", "derived", "crispulator_facs_repeated_metrics.csv"
  ),
  row.names = FALSE
)
write.csv(
  baseline_concordance,
  file.path(
    "data", "derived", "crispulator_facs_design_concordance.csv"
  ),
  row.names = FALSE
)
write.csv(
  baseline_summary,
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

# B: baseline cross-seed average precision and directional recall.
sorted <- baseline_metrics[baseline_metrics$method != "Bulk vs input", ]
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
bulk <- baseline_metrics[baseline_metrics$method == "Bulk vs input", ]
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
set.seed(20250724)
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

# Direct, seed-paired method differences for the primary three-sample design.
primary <- all_metrics[
  all_metrics$design == "Low + bulk + high" &
    all_metrics$method %in% c("BARCS", "MAGeCK-MLE"),
]
comparison_keys <- c(
  "scenario_id", "seed", "moi", "high_quality_guide_fraction",
  "genes", "is_baseline"
)
comparison_metrics <- c(
  "average_precision",
  "directional_recall_fdr_0_10",
  "f1_fdr_0_10"
)
barcs <- primary[
  primary$method == "BARCS",
  c(comparison_keys, comparison_metrics)
]
mageck <- primary[
  primary$method == "MAGeCK-MLE",
  c(comparison_keys, comparison_metrics)
]
method_differences <- merge(
  barcs,
  mageck,
  by = comparison_keys,
  suffixes = c("_barcs", "_mageck")
)
for (metric in comparison_metrics) {
  method_differences[[paste0(metric, "_difference")]] <-
    method_differences[[paste0(metric, "_barcs")]] -
    method_differences[[paste0(metric, "_mageck")]]
}
write.csv(
  method_differences,
  file.path(
    "data", "derived",
    "crispulator_facs_parameter_grid_method_differences.csv"
  ),
  row.names = FALSE
)

# Positive values favor BARCS and negative values favor MAGeCK-MLE. A full
# factorial plot marginalizes over the other two parameters; the default
# one-at-a-time plot holds them at their baseline values.
parameter_specs <- list(
  list(column = "moi", label = "MOI", panel = "A"),
  list(
    column = "high_quality_guide_fraction",
    label = "High-quality guide fraction",
    panel = "B"
  ),
  list(column = "genes", label = "Number of genes", panel = "C")
)
difference_specs <- list(
  list(
    column = "average_precision_difference",
    label = "Average precision",
    colour = "#0072B2",
    pch = 16,
    lty = 1
  ),
  list(
    column = "directional_recall_fdr_0_10_difference",
    label = "Directional recall",
    colour = "#D55E00",
    pch = 17,
    lty = 2
  ),
  list(
    column = "f1_fdr_0_10_difference",
    label = "F1",
    colour = "#009E73",
    pch = 15,
    lty = 3
  )
)

summarize_sweep <- function(parameter) {
  data <- method_differences
  if (grid_mode == "one_at_a_time") {
    fixed <- setdiff(
      c("moi", "high_quality_guide_fraction", "genes"),
      parameter
    )
    baseline_values <- c(
      moi = baseline_moi,
      high_quality_guide_fraction = baseline_quality,
      genes = baseline_genes
    )
    keep <- rep(TRUE, nrow(data))
    for (column in fixed) {
      keep <- keep & data[[column]] == baseline_values[[column]]
    }
    data <- data[keep, ]
  }
  do.call(rbind, lapply(
    split(data, data[[parameter]]),
    function(group) {
      row <- data.frame(value = group[[parameter]][1L])
      for (spec in difference_specs) {
        values <- group[[spec$column]]
        row[[paste0(spec$column, "_mean")]] <- mean(values)
        row[[paste0(spec$column, "_se")]] <-
          sd(values) / sqrt(length(values))
      }
      row
    }
  ))
}

sweep_summaries <- lapply(
  parameter_specs,
  function(spec) summarize_sweep(spec$column)
)
all_limits <- unlist(lapply(sweep_summaries, function(summary) {
  unlist(lapply(difference_specs, function(spec) {
    estimate <- summary[[paste0(spec$column, "_mean")]]
    error <- summary[[paste0(spec$column, "_se")]]
    c(estimate - error, estimate + error)
  }))
}))
y_limit <- range(c(-0.05, 0.05, all_limits), finite = TRUE)
y_padding <- max(0.01, diff(y_limit) * 0.08)
y_limit <- y_limit + c(-y_padding, y_padding)

pdf(
  file.path("figures", "crispulator_facs_parameter_sensitivity.pdf"),
  width = 9,
  height = 3.4,
  useDingbats = FALSE
)
par(mfrow = c(1, 3), mar = c(4.3, 4.3, 2.4, 0.8))
for (index in seq_along(parameter_specs)) {
  parameter <- parameter_specs[[index]]
  summary <- sweep_summaries[[index]]
  summary <- summary[order(summary$value), ]
  plot(
    range(summary$value),
    y_limit,
    type = "n",
    xlab = parameter$label,
    ylab = if (index == 1L) "BARCS minus MAGeCK-MLE" else "",
    main = paste(parameter$panel, parameter$label),
    bty = "l"
  )
  abline(h = 0, lty = 3, col = "#666666")
  for (spec in difference_specs) {
    estimate <- summary[[paste0(spec$column, "_mean")]]
    error <- summary[[paste0(spec$column, "_se")]]
    lines(
      summary$value,
      estimate,
      col = spec$colour,
      lwd = 1.6,
      lty = spec$lty
    )
    points(
      summary$value,
      estimate,
      col = spec$colour,
      bg = "white",
      pch = spec$pch
    )
    has_error <- is.finite(error) & error > 0
    arrows(
      summary$value[has_error],
      (estimate - error)[has_error],
      summary$value[has_error],
      (estimate + error)[has_error],
      angle = 90,
      code = 3,
      length = 0.035,
      col = spec$colour
    )
  }
  if (index == 2L) {
    legend(
      "bottom",
      legend = vapply(difference_specs, `[[`, character(1L), "label"),
      col = vapply(difference_specs, `[[`, character(1L), "colour"),
      pch = vapply(difference_specs, `[[`, numeric(1L), "pch"),
      lty = vapply(difference_specs, `[[`, numeric(1L), "lty"),
      bty = "n",
      cex = 0.72
    )
  }
}
dev.off()

print(summary_rows)
