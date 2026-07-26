#!/usr/bin/env Rscript

# Repeat the three-method BARCS CRISPulator benchmark across prespecified
# parameters and seeds. Per-run count matrices remain under git-ignored
# `results/`; compact metrics and figures are versioned.

options(stringsAsFactors = FALSE)
analysis_protocol <- "barcs-three-methods-v1"
expected_methods <- c(
  "BARCS-original", "BARCS-partial", "BARCS-EB"
)

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

cache_scenario_name <- function(moi, quality, genes, replicates) {
  name <- sprintf(
    "moi_%s_quality_%s",
    format_parameter(moi),
    format_parameter(quality)
  )
  if (genes != 400L) {
    name <- paste0(name, "_genes_", genes)
  }
  if (replicates != 4L) {
    name <- paste0(name, "_replicates_", replicates)
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
baseline_replicates <- as.integer(Sys.getenv("CRISPULATOR_REPLICATES", "4"))
moi_values <- parse_numeric_values(
  "CRISPULATOR_MOI_VALUES", c(0.10, 0.25, 0.40)
)
quality_values <- parse_numeric_values(
  "CRISPULATOR_HIGH_QUALITY_GUIDE_FRACTION_VALUES", c(0.60, 0.75, 0.90)
)
gene_values <- as.integer(parse_numeric_values(
  "CRISPULATOR_GENE_VALUES", c(100, 400, 1000)
))
replicate_values <- as.integer(parse_numeric_values(
  "CRISPULATOR_REPLICATE_VALUES", c(1, 3, 4, 6)
))
grid_mode <- Sys.getenv("CRISPULATOR_GRID_MODE", "one_at_a_time")

if (!is.finite(baseline_moi) ||
    baseline_moi <= 0 || baseline_moi >= 0.5 ||
    any(moi_values <= 0 | moi_values >= 0.5)) {
  stop("MOI values must be greater than zero and below 0.5.")
}
if (!is.finite(baseline_quality) ||
    baseline_quality < 0 || baseline_quality > 1 ||
    any(quality_values < 0 | quality_values > 1)) {
  stop("Guide-quality fractions must be between zero and one.")
}
if (!is.finite(baseline_genes) || baseline_genes < 20L ||
    any(gene_values < 20L)) {
  stop("Gene counts must be integers of at least 20.")
}
if (!is.finite(baseline_replicates) || baseline_replicates < 1L ||
    any(replicate_values < 1L)) {
  stop("Replicate counts must be positive integers.")
}
if (!grid_mode %in% c("single", "one_at_a_time", "full_factorial")) {
  stop(
    "`CRISPULATOR_GRID_MODE` must be single, one_at_a_time, or full_factorial."
  )
}

baseline <- data.frame(
  moi = baseline_moi,
  high_quality_guide_fraction = baseline_quality,
  genes = baseline_genes,
  replicates = baseline_replicates
)
scenarios <- switch(
  grid_mode,
  single = baseline,
  one_at_a_time = unique(rbind(
    baseline,
    data.frame(
      moi = moi_values,
      high_quality_guide_fraction = baseline_quality,
      genes = baseline_genes,
      replicates = baseline_replicates
    ),
    data.frame(
      moi = baseline_moi,
      high_quality_guide_fraction = quality_values,
      genes = baseline_genes,
      replicates = baseline_replicates
    ),
    data.frame(
      moi = baseline_moi,
      high_quality_guide_fraction = baseline_quality,
      genes = gene_values,
      replicates = baseline_replicates
    ),
    data.frame(
      moi = baseline_moi,
      high_quality_guide_fraction = baseline_quality,
      genes = baseline_genes,
      replicates = replicate_values
    )
  )),
  full_factorial = expand.grid(
    moi = moi_values,
    high_quality_guide_fraction = quality_values,
    genes = gene_values,
    replicates = replicate_values,
    KEEP.OUT.ATTRS = FALSE
  )
)
scenarios$scenario_id <- sprintf(
  "moi_%s_quality_%s_genes_%d_replicates_%d",
  vapply(scenarios$moi, format_parameter, character(1L)),
  vapply(
    scenarios$high_quality_guide_fraction,
    format_parameter,
    character(1L)
  ),
  scenarios$genes,
  scenarios$replicates
)
scenarios$is_baseline <- with(
  scenarios,
  moi == baseline_moi &
    high_quality_guide_fraction == baseline_quality &
    genes == baseline_genes &
    replicates == baseline_replicates
)
scenarios <- scenarios[
  order(
    !scenarios$is_baseline,
    scenarios$replicates,
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
    as.numeric(values[["multiplicity_of_infection"]]), run$moi
  )) &&
    isTRUE(all.equal(
      as.numeric(values[["high_quality_guide_fraction"]]),
      run$high_quality_guide_fraction
    )) &&
    identical(as.integer(values[["genes"]]), as.integer(run$genes)) &&
    identical(
      as.integer(values[["replicates"]]),
      as.integer(run$replicates)
    ) &&
    identical(as.integer(values[["seed"]]), as.integer(run$seed))
}

analysis_matches <- function(metric_file, replicates) {
  if (!file.exists(metric_file)) {
    return(FALSE)
  }
  metrics <- tryCatch(read.csv(metric_file), error = function(error) NULL)
  !is.null(metrics) &&
    identical(sort(unique(metrics$method)), sort(expected_methods)) &&
    all(metrics$analysis_protocol == analysis_protocol) &&
    (
      replicates >= 2L ||
        identical(unique(metrics$design), "Low + bulk + high")
    )
}

all_metrics <- lapply(seq_len(nrow(runs)), function(index) {
  run <- runs[index, ]
  cache_scenario <- cache_scenario_name(
    run$moi,
    run$high_quality_guide_fraction,
    run$genes,
    run$replicates
  )
  run_dir <- file.path(
    root, cache_scenario, paste0("seed_", run$seed)
  )
  data_dir <- file.path(run_dir, "data")
  result_dir <- file.path(run_dir, "analysis_three_methods")
  figure_path <- file.path(run_dir, "three_methods_benchmark.pdf")
  metric_file <- file.path(result_dir, "benchmark_metrics.csv")
  parameter_file <- file.path(data_dir, "parameters.tsv")
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  rerun <- identical(Sys.getenv("RERUN_CRISPULATOR"), "1")
  simulation_needed <- !cache_matches(parameter_file, run) || rerun
  analysis_needed <- simulation_needed ||
    !analysis_matches(metric_file, run$replicates) || rerun
  if (simulation_needed || analysis_needed) {
    message(sprintf(
      paste0(
        "[%d/%d] MOI=%s quality=%s genes=%d ",
        "replicates=%d seed=%d"
      ),
      index,
      nrow(runs),
      format_parameter(run$moi),
      format_parameter(run$high_quality_guide_fraction),
      run$genes,
      run$replicates,
      run$seed
    ))
  }
  if (simulation_needed) {
    status <- system2(
      "julia",
      c(
        "--project=julia",
        file.path("julia", "simulate_crispulator_facs.jl"),
        data_dir,
        as.character(run$replicates),
        as.character(run$seed),
        as.character(run$moi),
        as.character(run$high_quality_guide_fraction),
        as.character(run$genes)
      )
    )
    if (status != 0) {
      stop("CRISPulator simulation failed for ", run$scenario_id, ".")
    }
  }
  if (analysis_needed) {
    status <- system2(
      "Rscript",
      file.path("examples", "crispulator_facs_benchmark.R"),
      env = c(
        paste0("CRISPULATOR_DATA_DIR=", data_dir),
        paste0("CRISPULATOR_RESULT_DIR=", result_dir),
        paste0("CRISPULATOR_FIGURE_PATH=", figure_path),
        paste0("BARCS_NCORES=", Sys.getenv("BARCS_NCORES", "4"))
      )
    )
    if (status != 0) {
      stop("Three-method BARCS analysis failed for ", run$scenario_id, ".")
    }
  }
  metrics <- read.csv(metric_file)
  metrics$scenario_id <- run$scenario_id
  metrics$seed <- run$seed
  metrics$moi <- run$moi
  metrics$high_quality_guide_fraction <-
    run$high_quality_guide_fraction
  metrics$genes <- run$genes
  metrics$replicates <- run$replicates
  metrics$is_baseline <- run$is_baseline
  metrics
})
all_metrics <- do.call(rbind, all_metrics)
write.csv(
  all_metrics,
  file.path(
    "data", "derived", "crispulator_facs_repeated_metrics.csv"
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

runtime_rows <- do.call(rbind, lapply(seq_len(nrow(runs)), function(index) {
  run <- runs[index, ]
  cache_scenario <- cache_scenario_name(
    run$moi,
    run$high_quality_guide_fraction,
    run$genes,
    run$replicates
  )
  runtime_path <- file.path(
    root,
    cache_scenario,
    paste0("seed_", run$seed),
    "analysis_three_methods",
    "low_bulk_high_runtime.csv"
  )
  runtime <- read.csv(runtime_path)
  runtime$scenario_id <- run$scenario_id
  runtime$seed <- run$seed
  runtime$moi <- run$moi
  runtime$high_quality_guide_fraction <-
    run$high_quality_guide_fraction
  runtime$genes <- run$genes
  runtime$replicates <- run$replicates
  runtime$total_seconds <-
    runtime$shared_guide_fit_seconds +
    runtime$gene_aggregation_seconds
  runtime
}))
write.csv(
  runtime_rows,
  file.path(
    "data", "derived", "crispulator_facs_multimethod_runtime.csv"
  ),
  row.names = FALSE
)

summarize_metric <- function(data, grouping, metric) {
  groups <- split(
    seq_len(nrow(data)),
    interaction(data[grouping], drop = TRUE, lex.order = TRUE)
  )
  do.call(rbind, lapply(groups, function(index) {
    group <- data[index, , drop = FALSE]
    values <- group[[metric]]
    row <- group[1L, grouping, drop = FALSE]
    row$metric <- metric
    row$n <- sum(is.finite(values))
    row$mean <- mean(values, na.rm = TRUE)
    row$sd <- sd(values, na.rm = TRUE)
    row$se <- row$sd / sqrt(row$n)
    row
  }))
}

summary_metrics <- c(
  "auroc",
  "average_precision",
  "effect_spearman_active",
  "direction_accuracy_active",
  "directional_recall_fdr_0_10",
  "empirical_fdp_fdr_0_10",
  "f1_fdr_0_10",
  "negative_control_p_below_0_05"
)
summary_grouping <- c(
  "scenario_id", "method", "design", "moi",
  "high_quality_guide_fraction", "genes", "replicates",
  "is_baseline"
)
parameter_summary <- do.call(rbind, lapply(summary_metrics, function(metric) {
  summarize_metric(all_metrics, summary_grouping, metric)
}))
write.csv(
  parameter_summary,
  file.path(
    "data", "derived", "crispulator_facs_parameter_grid_summary.csv"
  ),
  row.names = FALSE
)
write.csv(
  parameter_summary,
  file.path(
    "data", "derived", "crispulator_facs_repeated_summary.csv"
  ),
  row.names = FALSE
)

runtime_summary <- summarize_metric(
  runtime_rows,
  c("scenario_id", "method", "moi", "high_quality_guide_fraction",
    "genes", "replicates"),
  "total_seconds"
)
write.csv(
  runtime_summary,
  file.path(
    "data", "derived",
    "crispulator_facs_multimethod_runtime_summary.csv"
  ),
  row.names = FALSE
)

primary <- all_metrics[
  all_metrics$design == "Low + bulk + high",
  ,
  drop = FALSE
]
method_pairs <- list(
  c("BARCS-partial", "BARCS-original"),
  c("BARCS-EB", "BARCS-original"),
  c("BARCS-EB", "BARCS-partial")
)
paired_rows <- list()
pair_index <- 0L
for (metric in summary_metrics) {
  for (pair in method_pairs) {
    left <- primary[
      primary$method == pair[1L],
      c("scenario_id", "seed", metric)
    ]
    right <- primary[
      primary$method == pair[2L],
      c("scenario_id", "seed", metric)
    ]
    names(left)[3L] <- "left"
    names(right)[3L] <- "right"
    paired <- merge(left, right, by = c("scenario_id", "seed"))
    paired$difference <- paired$left - paired$right
    groups <- split(paired, paired$scenario_id)
    for (scenario in names(groups)) {
      group <- groups[[scenario]]
      pair_index <- pair_index + 1L
      paired_rows[[pair_index]] <- data.frame(
        scenario_id = scenario,
        metric = metric,
        method = pair[1L],
        reference = pair[2L],
        n = nrow(group),
        mean_difference = mean(group$difference),
        sd_difference = sd(group$difference),
        lower_95 = mean(group$difference) -
          qt(0.975, df = nrow(group) - 1L) *
            sd(group$difference) / sqrt(nrow(group)),
        upper_95 = mean(group$difference) +
          qt(0.975, df = nrow(group) - 1L) *
            sd(group$difference) / sqrt(nrow(group))
      )
    }
  }
}
paired_summary <- do.call(rbind, paired_rows)
write.csv(
  paired_summary,
  file.path(
    "data", "derived",
    "crispulator_facs_three_method_paired_differences.csv"
  ),
  row.names = FALSE
)
write.csv(
  paired_summary,
  file.path(
    "data", "derived",
    "crispulator_facs_parameter_grid_method_differences.csv"
  ),
  row.names = FALSE
)

replicate_metrics <- primary[
  primary$moi == baseline_moi &
    primary$high_quality_guide_fraction == baseline_quality &
    primary$genes == baseline_genes,
  ,
  drop = FALSE
]
replicate_summary <- do.call(rbind, lapply(
  c("average_precision", "directional_recall_fdr_0_10",
    "empirical_fdp_fdr_0_10"),
  function(metric) {
    summarize_metric(
      replicate_metrics,
      c("method", "replicates"),
      metric
    )
  }
))
write.csv(
  replicate_summary,
  file.path(
    "data", "derived",
    "crispulator_facs_multimethod_replicate_summary.csv"
  ),
  row.names = FALSE
)

method_colours <- c(
  `BARCS-original` = "#0072B2",
  `BARCS-partial` = "#009E73",
  `BARCS-EB` = "#D55E00"
)
replicate_figure <- file.path(
  "figures", "crispulator_facs_multimethod_replicates.pdf"
)
pdf(replicate_figure, width = 10.5, height = 3.8, useDingbats = FALSE)
layout(
  matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE),
  heights = c(1, 0.14)
)
figure_metrics <- c(
  average_precision = "Average precision",
  directional_recall_fdr_0_10 = "Directional recall",
  empirical_fdp_fdr_0_10 = "Realized FDP"
)
for (metric in names(figure_metrics)) {
  par(mar = c(4.2, 4.2, 2.5, 0.8))
  panel <- replicate_summary[
    replicate_summary$metric == metric, , drop = FALSE
  ]
  ylim <- range(c(0, panel$mean + panel$se), finite = TRUE)
  if (metric != "empirical_fdp_fdr_0_10") {
    ylim <- c(0, 1)
  }
  plot(
    range(panel$replicates),
    ylim,
    type = "n",
    xlab = "Screen replicates",
    ylab = figure_metrics[[metric]],
    main = figure_metrics[[metric]],
    bty = "l"
  )
  for (method in expected_methods) {
    rows <- panel[panel$method == method, ]
    rows <- rows[order(rows$replicates), ]
    lines(
      rows$replicates,
      rows$mean,
      type = "b",
      pch = 16,
      lwd = 2,
      col = method_colours[[method]]
    )
  }
  if (metric == "empirical_fdp_fdr_0_10") {
    abline(h = 0.10, lty = 2, col = "#555555")
  }
}
par(mar = c(0, 0, 0, 0))
plot.new()
legend(
  "center",
  legend = expected_methods,
  col = unname(method_colours[expected_methods]),
  lwd = 2,
  pch = 16,
  horiz = TRUE,
  bty = "n"
)
dev.off()

print(parameter_summary)
print(paired_summary)
