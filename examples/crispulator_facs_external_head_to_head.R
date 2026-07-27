#!/usr/bin/env Rscript

# Head-to-head comparison of the four BARCS gene statistics with four
# established count-analysis methods. External fits are read from the
# previously rerun `analysis/` directories only after the historical BARCS
# effect vector is verified against the current BARCS-original result for
# every run. All metrics are then recomputed with one current evaluator.
#
# This script writes the per-run metric and provenance tables for every
# scenario and then sources the aggregation stage, which restricts the
# headline summary to the multi-replicate settings and reports the
# one-replicate boundary separately. See
# `examples/crispulator_facs_external_head_to_head_aggregate.R`.

options(stringsAsFactors = FALSE)
analysis_protocol <- "barcs-external-headtohead-v1"
seeds <- 20250724L:20250728L
scenario_table <- read.csv(
  file.path("data", "derived", "crispulator_facs_parameter_grid.csv")
)
expected_methods <- c(
  "BARCS-original", "BARCS-NORM", "BARCS-partial", "BARCS-EB",
  "MAGeCK-MLE", "edgeR-QL", "DESeq2", "limma-voom"
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

auroc <- function(truth, score) {
  score_rank <- rank(score, ties.method = "average")
  n_positive <- sum(truth)
  n_negative <- sum(!truth)
  (
    sum(score_rank[truth]) -
      n_positive * (n_positive + 1) / 2
  ) / (n_positive * n_negative)
}

average_precision <- function(truth, score) {
  ordered_truth <- truth[order(score, decreasing = TRUE)]
  precision_at_rank <- cumsum(ordered_truth) / seq_along(ordered_truth)
  mean(precision_at_rank[ordered_truth])
}

evaluate <- function(method, gene_result, gene_truth) {
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
  score <- -log10(pmax(assessed$p_value, .Machine$double.xmin))
  called <- assessed$fdr < 0.10
  active <- assessed$active
  sign_match <- sign(assessed$estimate) == assessed$expected_sign
  true_positive <- sum(called & active)
  false_positive <- sum(called & !active)
  false_negative <- sum(!called & active)
  precision <- if ((true_positive + false_positive) > 0) {
    true_positive / (true_positive + false_positive)
  } else {
    0
  }
  recall <- if ((true_positive + false_negative) > 0) {
    true_positive / (true_positive + false_negative)
  } else {
    0
  }
  f1 <- if ((precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    0
  }
  negative_control <- assessed$class == "negcontrol"
  data.frame(
    analysis_protocol = analysis_protocol,
    method = method,
    genes = nrow(assessed),
    active_genes = sum(active),
    auroc = auroc(active, score),
    average_precision = average_precision(active, score),
    effect_spearman_active = cor(
      assessed$estimate[active],
      assessed$theoretical_phenotype[active],
      method = "spearman"
    ),
    direction_accuracy_active = mean(sign_match[active]),
    calls_fdr_0_10 = sum(called),
    directional_recall_fdr_0_10 =
      mean(called[active] & sign_match[active]),
    empirical_fdp_fdr_0_10 =
      if (sum(called) > 0) false_positive / sum(called) else 0,
    f1_fdr_0_10 = f1,
    negative_control_p_below_0_05 =
      mean(assessed$p_value[negative_control] < 0.05)
  )
}

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

metric_rows <- list()
provenance_rows <- list()
row_index <- 0L
provenance_index <- 0L
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
    data_dir <- file.path(run_root, "data")
    historical_dir <- file.path(run_root, "analysis")
    current_dir <- file.path(run_root, "analysis_three_methods")
    truth_path <- file.path(data_dir, "gene_truth.tsv")
    count_path <- file.path(data_dir, "counts.tsv")
    design_path <- file.path(data_dir, "sample_design.tsv")
    historical_barcs_path <- file.path(
      historical_dir, "barcs_low_bulk_high_gene_results.csv"
    )
    current_barcs_path <- file.path(
      current_dir, method_file[["BARCS-original"]]
    )
    required <- c(
      truth_path, count_path, design_path,
      historical_barcs_path,
      file.path(
        current_dir,
        method_file[c(
          "BARCS-original", "BARCS-NORM", "BARCS-partial", "BARCS-EB"
        )]
      ),
      file.path(
        historical_dir,
        method_file[c("MAGeCK-MLE", "edgeR-QL", "DESeq2", "limma-voom")]
      )
    )
    if (!all(file.exists(required))) {
      stop(
        "Missing head-to-head result files for ",
        scenario$scenario_id, ", seed ", seed, "."
      )
    }

    historical_barcs <- read.csv(historical_barcs_path)
    current_barcs <- read.csv(current_barcs_path)
    identity_check <- merge(
      historical_barcs[
        , c("gene", "estimate", "p_value", "fdr")
      ],
      current_barcs[, c("gene", "estimate", "p_value", "fdr")],
      by = "gene",
      suffixes = c("_historical", "_current")
    )
    if (nrow(identity_check) != nrow(current_barcs) ||
        nrow(identity_check) != nrow(historical_barcs)) {
      stop("Historical and current BARCS gene universes differ.")
    }
    max_effect_difference <- max(abs(
      identity_check$estimate_historical -
        identity_check$estimate_current
    ))
    max_p_difference <- max(abs(
      identity_check$p_value_historical -
        identity_check$p_value_current
    ))
    max_fdr_difference <- max(abs(
      identity_check$fdr_historical -
        identity_check$fdr_current
    ))
    if (max_effect_difference > 1e-12 ||
        max_p_difference > 1e-3 ||
        max_fdr_difference > 1e-3) {
      stop(
        "Historical external fits fail the BARCS provenance check for ",
        scenario$scenario_id, ", seed ", seed, "."
      )
    }

    gene_truth <- read.delim(truth_path)
    gene_truth$active <-
      tolower(as.character(gene_truth$active)) == "true"
    method_results <- list(
      `BARCS-original` = current_barcs,
      `BARCS-NORM` = read.csv(file.path(
        current_dir, method_file[["BARCS-NORM"]]
      )),
      `BARCS-partial` = read.csv(file.path(
        current_dir, method_file[["BARCS-partial"]]
      )),
      `BARCS-EB` = read.csv(file.path(
        current_dir, method_file[["BARCS-EB"]]
      )),
      `MAGeCK-MLE` = read.csv(file.path(
        historical_dir, method_file[["MAGeCK-MLE"]]
      )),
      `edgeR-QL` = read.csv(file.path(
        historical_dir, method_file[["edgeR-QL"]]
      )),
      DESeq2 = read.csv(file.path(
        historical_dir, method_file[["DESeq2"]]
      )),
      `limma-voom` = read.csv(file.path(
        historical_dir, method_file[["limma-voom"]]
      ))
    )
    source_paths <- c(
      file.path(
        current_dir,
        method_file[c(
          "BARCS-original", "BARCS-NORM", "BARCS-partial", "BARCS-EB"
        )]
      ),
      file.path(
        historical_dir,
        method_file[c("MAGeCK-MLE", "edgeR-QL", "DESeq2", "limma-voom")]
      )
    )
    names(source_paths) <- expected_methods
    for (method in expected_methods) {
      row_index <- row_index + 1L
      metric <- evaluate(method, method_results[[method]], gene_truth)
      metric$scenario_id <- scenario$scenario_id
      metric$seed <- seed
      metric$moi <- scenario$moi
      metric$high_quality_guide_fraction <-
        scenario$high_quality_guide_fraction
      metric$genes_simulated <- scenario$genes
      metric$replicates <- scenario$replicates
      metric$is_baseline <- scenario$is_baseline
      metric_rows[[row_index]] <- metric

      provenance_index <- provenance_index + 1L
      provenance_rows[[provenance_index]] <- data.frame(
        analysis_protocol = analysis_protocol,
        scenario_id = scenario$scenario_id,
        seed = seed,
        method = method,
        result_source = if (grepl("^BARCS", method)) {
          "current four-method run"
        } else {
          "verified stored external fit"
        },
        result_path = source_paths[[method]],
        result_md5 = unname(tools::md5sum(source_paths[[method]])),
        count_path = count_path,
        count_md5 = unname(tools::md5sum(count_path)),
        truth_path = truth_path,
        truth_md5 = unname(tools::md5sum(truth_path)),
        design_path = design_path,
        design_md5 = unname(tools::md5sum(design_path)),
        barcs_max_effect_difference = max_effect_difference,
        barcs_max_p_difference = max_p_difference,
        barcs_max_fdr_difference = max_fdr_difference
      )
    }
  }
}
metrics <- do.call(rbind, metric_rows)
provenance <- do.call(rbind, provenance_rows)
stopifnot(
  nrow(metrics) == nrow(scenario_table) * length(seeds) *
    length(expected_methods),
  identical(sort(unique(metrics$method)), sort(expected_methods)),
  all(metrics$analysis_protocol == analysis_protocol)
)
write.csv(
  metrics,
  file.path(
    "data", "derived",
    "crispulator_facs_external_head_to_head_metrics.csv"
  ),
  row.names = FALSE
)
write.csv(
  provenance,
  file.path(
    "data", "derived",
    "crispulator_facs_external_head_to_head_provenance.csv"
  ),
  row.names = FALSE
)

source(file.path(
  "examples", "crispulator_facs_external_head_to_head_aggregate.R"
))
