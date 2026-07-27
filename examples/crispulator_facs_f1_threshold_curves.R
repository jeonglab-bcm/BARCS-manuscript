#!/usr/bin/env Rscript

# F1 and calibration curves across nominal gene-FDR thresholds for the eight
# methods in the CRISPulator low + bulk + high comparison.
#
# Reporting F1 at a single cutoff invites a method to win by being
# miscalibrated. A method whose realized error rate runs well above what it
# advertises simply calls more genes, and at a fixed nominal FDR that buys
# recall and therefore F1 -- without the method being any better at
# distinguishing active genes from inactive ones.
#
# So the comparison is run across a range of nominal thresholds and two things
# are reported side by side: F1, and the realized false-discovery proportion
# against the requested one. Reading them together is what separates "detects
# more" from "asks for less". It is also what shows the F1 ordering here is
# fragile -- the leader changes with the cutoff -- which is stated in the
# manuscript rather than resolved in a method's favour.
#
# This is stage one of two. It fits and scores every method on every run and
# writes per-run rows; `crispulator_facs_f1_threshold_curves_aggregate.R`
# applies the scope rule and builds the summaries and the figure. The split
# means a change to what counts as in scope is a cheap rerun of stage two
# rather than a full refit.
#
# Reads the per-run screens written by the parameter-grid benchmark, so that
# has to have run first.
#
#     pixi run bench-facs-grid
#     pixi run bench-f1-curves
#     pixi run bench-f1-curves-aggregate
#
# The genome-scale version of this scan, over thirteen log-spaced cutoffs from
# 1e-6 to 0.20, lives in `examples/crispulator_facs_moi_10k_benchmark.R`.

options(stringsAsFactors = FALSE)
analysis_protocol <- "crispulator-f1-threshold-curves-v1"
thresholds <- c(0.10, 0.05, 0.01, 0.005, 0.001)
seeds <- 20250724L:20250728L
scenario_table <- read.csv(
  file.path("data", "derived", "crispulator_facs_parameter_grid.csv")
)
methods <- c(
  "BARCS-original", "BARCS-NORM", "BARCS-partial", "BARCS-EB",
  "MAGeCK-MLE", "edgeR-QL", "DESeq2", "limma-voom"
)
method_file <- c(
  `BARCS-original` = "low_bulk_high_original_gene_results.csv",
  `BARCS-NORM` = "low_bulk_high_normal_gene_results.csv",
  `BARCS-partial` = "low_bulk_high_partial_gene_results.csv",
  `BARCS-EB` = "low_bulk_high_eb_gene_results.csv",
  `MAGeCK-MLE` = "mageck_mle_low_bulk_high_gene_results.csv",
  `edgeR-QL` = "edger_ql_low_bulk_high_gene_results.csv",
  DESeq2 = "deseq2_low_bulk_high_gene_results.csv",
  `limma-voom` = "limma_voom_low_bulk_high_gene_results.csv"
)

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

safe_ratio <- function(numerator, denominator, empty = NA_real_) {
  if (denominator > 0) numerator / denominator else empty
}

evaluate_threshold <- function(
  method,
  gene_result,
  gene_truth,
  threshold
) {
  assessed <- merge(
    gene_truth,
    gene_result[, c("gene", "estimate", "p_value", "fdr")],
    by = "gene",
    all.x = TRUE
  )
  assessed <- assessed[
    is.finite(assessed$estimate) &
      is.finite(assessed$p_value) &
      is.finite(assessed$fdr),
    ,
    drop = FALSE
  ]
  active <- assessed$active
  called <- assessed$fdr < threshold
  true_positive <- sum(called & active)
  false_positive <- sum(called & !active)
  true_negative <- sum(!called & !active)
  false_negative <- sum(!called & active)
  precision <- safe_ratio(
    true_positive,
    true_positive + false_positive,
    empty = 0
  )
  recall <- safe_ratio(
    true_positive,
    true_positive + false_negative,
    empty = 0
  )
  f1 <- safe_ratio(
    2 * precision * recall,
    precision + recall,
    empty = 0
  )
  sign_match <- sign(assessed$estimate) == assessed$expected_sign
  data.frame(
    analysis_protocol = analysis_protocol,
    method = method,
    nominal_fdr = threshold,
    genes = nrow(assessed),
    active_genes = sum(active),
    calls = sum(called),
    true_positive = true_positive,
    false_positive = false_positive,
    true_negative = true_negative,
    false_negative = false_negative,
    precision = precision,
    recall = recall,
    directional_recall = mean(called[active] & sign_match[active]),
    f1 = f1,
    realized_fdp = safe_ratio(
      false_positive,
      true_positive + false_positive,
      empty = 0
    ),
    specificity = safe_ratio(
      true_negative,
      true_negative + false_positive
    ),
    balanced_accuracy = mean(c(
      recall,
      safe_ratio(true_negative, true_negative + false_positive)
    ))
  )
}

rows <- list()
row_index <- 0L
for (scenario_index in seq_len(nrow(scenario_table))) {
  scenario <- scenario_table[scenario_index, ]
  cache_scenario <- cache_scenario_name(
    scenario$moi,
    scenario$high_quality_guide_fraction,
    scenario$genes,
    scenario$replicates
  )
  for (seed in seeds) {
    run_root <- file.path(
      "results", "crispulator_facs", "repeated",
      cache_scenario, paste0("seed_", seed)
    )
    truth_path <- file.path(run_root, "data", "gene_truth.tsv")
    current_dir <- file.path(run_root, "analysis_three_methods")
    external_dir <- file.path(run_root, "analysis")
    result_paths <- c(
      file.path(
        current_dir,
        method_file[c(
          "BARCS-original", "BARCS-NORM", "BARCS-partial", "BARCS-EB"
        )]
      ),
      file.path(
        external_dir,
        method_file[c(
          "MAGeCK-MLE", "edgeR-QL", "DESeq2", "limma-voom"
        )]
      )
    )
    names(result_paths) <- methods
    required <- c(truth_path, result_paths)
    if (!all(file.exists(required))) {
      stop(
        "Missing truth or gene-result file(s) for ",
        scenario$scenario_id, ", seed ", seed, "."
      )
    }

    gene_truth <- read.delim(truth_path)
    gene_truth$active <-
      tolower(as.character(gene_truth$active)) == "true"
    method_results <- lapply(result_paths, read.csv)
    for (method in methods) {
      for (threshold in thresholds) {
        row_index <- row_index + 1L
        metric <- evaluate_threshold(
          method,
          method_results[[method]],
          gene_truth,
          threshold
        )
        metric$scenario_id <- scenario$scenario_id
        metric$seed <- seed
        metric$moi <- scenario$moi
        metric$high_quality_guide_fraction <-
          scenario$high_quality_guide_fraction
        metric$genes_simulated <- scenario$genes
        metric$replicates <- scenario$replicates
        metric$is_baseline <- scenario$is_baseline
        rows[[row_index]] <- metric
      }
    }
  }
}
metrics <- do.call(rbind, rows)
stopifnot(
  nrow(metrics) ==
    nrow(scenario_table) * length(seeds) * length(methods) *
      length(thresholds),
  setequal(unique(metrics$method), methods),
  setequal(unique(metrics$nominal_fdr), thresholds),
  all(metrics$f1 >= 0 & metrics$f1 <= 1),
  all(metrics$realized_fdp >= 0 & metrics$realized_fdp <= 1)
)
write.csv(
  metrics,
  file.path(
    "data", "derived", "crispulator_facs_f1_by_fdr.csv"
  ),
  row.names = FALSE
)

source(file.path(
  "examples", "crispulator_facs_f1_threshold_curves_aggregate.R"
))
