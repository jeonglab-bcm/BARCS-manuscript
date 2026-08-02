#!/usr/bin/env Rscript

# Held-out validation of the BARCS guide-dispersion moderation.
#
# The moderation defaults and the `qq_slope` null-scale band were chosen on a
# development set: the baseline scenario at the five seeds already used
# throughout this repository. Reporting on that same set would measure the
# choice, not the method. This script therefore regenerates the nine supported
# multi-replicate scenarios at three seeds that were never inspected during
# development, and evaluates every method once on them.
#
# Stage one needs Julia and the pinned environment in `julia/`:
#
#     Rscript examples/crispulator_facs_improved_barcs_holdout.R --simulate
#
# Stage two fits and evaluates, and needs edgeR, limma, and DESeq2:
#
#     Rscript examples/crispulator_facs_improved_barcs_holdout.R
#
# Simulations land under the ignored `results/` tree; only the metric tables
# are committed.

options(stringsAsFactors = FALSE)
analysis_protocol <- "crispulator-improved-barcs-holdout-v1"
nominal_fdr <- 0.10
# Two independent fresh-seed splits. The held-out split decided one default
# (see docs/barcs-external-method-comparison.md); the confirmatory split was
# generated and run afterwards to test that decision on data it did not touch.
seed_splits <- list(
  `held-out` = c(20260101L, 20260102L, 20260103L),
  confirmatory = c(20260201L, 20260202L, 20260203L)
)
simulation_root <- file.path(
  "results", "crispulator_facs", "improved_holdout"
)

scenarios <- data.frame(
  replicates = c(4L, 3L, 6L, 4L, 4L, 4L, 4L, 4L, 4L),
  moi = c(0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.10, 0.40, 0.25),
  high_quality_guide_fraction =
    c(0.9, 0.9, 0.9, 0.9, 0.6, 0.75, 0.9, 0.9, 0.9),
  genes = c(400L, 400L, 400L, 100L, 400L, 400L, 400L, 400L, 1000L)
)
# `format()` pads a vector to a common width, so format each value on its own.
compact_number <- function(value) {
  vapply(value, function(x) format(x, trim = TRUE, scientific = FALSE),
         character(1))
}
scenarios$scenario_id <- sprintf(
  "moi%s_q%s_g%d_r%d",
  compact_number(scenarios$moi),
  compact_number(scenarios$high_quality_guide_fraction),
  scenarios$genes, scenarios$replicates
)

run_directory <- function(scenario_id, seed) {
  file.path(simulation_root, sprintf("%s__seed_%d", scenario_id, seed))
}

if ("--simulate" %in% commandArgs(trailingOnly = TRUE)) {
  for (index in seq_len(nrow(scenarios))) {
    scenario <- scenarios[index, ]
    for (seed in unlist(seed_splits, use.names = FALSE)) {
      target <- run_directory(scenario$scenario_id, seed)
      if (file.exists(file.path(target, "counts.tsv"))) {
        next
      }
      status <- system2("julia", c(
        "--project=julia", file.path("julia", "simulate_crispulator_facs.jl"),
        target, scenario$replicates, seed, scenario$moi,
        scenario$high_quality_guide_fraction, scenario$genes
      ))
      if (status != 0L) {
        stop("Simulation failed for ", target, call. = FALSE)
      }
      message("simulated ", basename(target))
    }
  }
  message("Simulation complete. Re-run without --simulate to evaluate.")
  quit(save = "no")
}

needed <- c("edgeR", "limma", "DESeq2")
absent <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(absent)) {
  stop(
    "This comparison needs ", paste(absent, collapse = ", "),
    ". Install them from Bioconductor, or on Debian/Ubuntu with ",
    "apt-get install r-bioc-edger r-bioc-limma r-bioc-deseq2.",
    call. = FALSE
  )
}
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(DESeq2)
})
source(file.path("R", "load_barcs.R"))

auroc <- function(truth, score) {
  ranks <- rank(score, ties.method = "average")
  n_positive <- sum(truth)
  (sum(ranks[truth]) - n_positive * (n_positive + 1) / 2) /
    (n_positive * sum(!truth))
}
average_precision <- function(truth, score) {
  ordered_truth <- truth[order(score, decreasing = TRUE)]
  mean((cumsum(ordered_truth) / seq_along(ordered_truth))[ordered_truth])
}

fit_one_run <- function(directory) {
  count_table <- read.delim(
    file.path(directory, "counts.tsv"), check.names = FALSE
  )
  sample_data <- read.delim(
    file.path(directory, "sample_design.tsv"), check.names = FALSE
  )
  guide_truth <- read.delim(
    file.path(directory, "guide_truth.tsv"), check.names = FALSE
  )
  gene_truth <- read.delim(
    file.path(directory, "gene_truth.tsv"), check.names = FALSE
  )
  gene_truth$active <- tolower(as.character(gene_truth$active)) == "true"
  counts <- as.matrix(count_table[, sample_data$sample, drop = FALSE])
  storage.mode(counts) <- "double"

  conditional_normal_mean <- function(lower, upper) {
    mass <- pnorm(upper) - pnorm(lower)
    (dnorm(lower) - dnorm(upper)) / mass
  }
  q25 <- qnorm(0.25)
  sample_data$phenotype_z <- NA_real_
  sample_data$phenotype_z[sample_data$sample_type == "low"] <-
    conditional_normal_mean(-Inf, q25)
  sample_data$phenotype_z[sample_data$sample_type == "high"] <-
    conditional_normal_mean(-q25, Inf)
  sample_data$phenotype_z[sample_data$sample_type == "bulk"] <- 0
  sample_data$replicate <- factor(sample_data$replicate)

  keep_sample <- sample_data$sample_type %in% c("low", "bulk", "high")
  design_data <- droplevels(sample_data[keep_sample, , drop = FALSE])
  y <- counts[, keep_sample, drop = FALSE]
  formula <- ~ phenotype_z + replicate
  model_matrix <- model.matrix(formula, data = design_data)
  keep_guide <- rowSums(y) >= 30

  guide_result <- bb_screen(
    counts = y, totals = colSums(y), data = design_data,
    formula = formula, term = "phenotype_z",
    guide = guide_truth$guide, gene = guide_truth$gene,
    min_total_count = 30,
    ncores = as.integer(Sys.getenv("BARCS_NCORES", "4"))
  )
  negative_control <- guide_truth$class == "negcontrol"
  # Standard two-sided moderation, and the conservative one-way variant that
  # neither lowers a guide variance nor claims the prior degrees of freedom.
  moderated <- bb_moderate_dispersion(
    guide_result, trend = TRUE, one_way = FALSE, borrow_df = TRUE
  )
  conservative <- bb_moderate_dispersion(
    guide_result, trend = TRUE, one_way = TRUE, borrow_df = FALSE
  )

  gene_results <- list(
    `BARCS-original` = bb_gene_original(
      bb_calibrate_controls(guide_result, negative_control, alpha = 0.05),
      min_guides = 1L
    ),
    `BARCS-moderated` = bb_gene_original(
      bb_calibrate_controls(
        moderated, negative_control, alpha = 0.05, method = "qq_slope"
      ),
      min_guides = 1L
    ),
    `BARCS-moderated-one-way` = bb_gene_original(
      bb_calibrate_controls(
        conservative, negative_control, alpha = 0.05, method = "qq_slope"
      ),
      min_guides = 1L
    ),
    `BARCS-moderation-only` = bb_gene_original(
      bb_calibrate_controls(moderated, negative_control, alpha = 0.05),
      min_guides = 1L
    ),
    `BARCS-qq-slope-only` = bb_gene_original(
      bb_calibrate_controls(
        guide_result, negative_control, alpha = 0.05, method = "qq_slope"
      ),
      min_guides = 1L
    )
  )

  dge <- calcNormFactors(DGEList(counts = y[keep_guide, , drop = FALSE]))
  dge <- estimateDisp(dge, model_matrix, robust = TRUE)
  quasi_likelihood <- glmQLFTest(
    glmQLFit(dge, model_matrix, robust = TRUE), coef = "phenotype_z"
  )
  gene_results[["edgeR-QL"]] <- bb_gene_original(data.frame(
    gene = guide_truth$gene[keep_guide],
    estimate = quasi_likelihood$table$logFC,
    p_value = quasi_likelihood$table$PValue
  ), min_guides = 1L)

  voom_fit <- eBayes(lmFit(voom(dge, model_matrix), model_matrix))
  gene_results[["limma-voom"]] <- bb_gene_original(data.frame(
    gene = guide_truth$gene[keep_guide],
    estimate = voom_fit$coefficients[, "phenotype_z"],
    p_value = voom_fit$p.value[, "phenotype_z"]
  ), min_guides = 1L)

  deseq_result <- results(DESeq(DESeqDataSetFromMatrix(
    countData = round(y[keep_guide, , drop = FALSE]),
    colData = design_data, design = formula
  ), quiet = TRUE), name = "phenotype_z")
  gene_results[["DESeq2"]] <- bb_gene_original(data.frame(
    gene = guide_truth$gene[keep_guide],
    estimate = deseq_result$log2FoldChange,
    p_value = deseq_result$pvalue
  ), min_guides = 1L)

  do.call(rbind, lapply(names(gene_results), function(method) {
    assessed <- merge(
      gene_truth,
      gene_results[[method]][, c("gene", "estimate", "p_value", "fdr")],
      by = "gene", all.x = TRUE
    )
    assessed <- assessed[
      is.finite(assessed$estimate) & is.finite(assessed$p_value) &
        is.finite(assessed$fdr), , drop = FALSE
    ]
    score <- -log10(pmax(assessed$p_value, .Machine$double.xmin))
    called <- assessed$fdr < nominal_fdr
    active <- assessed$active
    sign_match <- sign(assessed$estimate) == assessed$expected_sign
    true_positive <- sum(called & active)
    false_positive <- sum(called & !active)
    false_negative <- sum(!called & active)
    data.frame(
      analysis_protocol = analysis_protocol,
      method = method,
      genes = nrow(assessed),
      active_genes = sum(active),
      auroc = auroc(active, score),
      average_precision = average_precision(active, score),
      calls = sum(called),
      directional_recall = mean(called[active] & sign_match[active]),
      realized_fdp = if (sum(called) > 0) {
        false_positive / sum(called)
      } else {
        0
      },
      f1 = if (true_positive > 0) {
        2 * true_positive /
          (2 * true_positive + false_positive + false_negative)
      } else {
        0
      },
      inactive_gene_p_below_0_05 = mean(assessed$p_value[!active] < 0.05),
      estimated_prior_df = attr(moderated, "prior_df")
    )
  }))
}

metric_rows <- list()
skipped <- character(0)
for (split in names(seed_splits)) {
  for (index in seq_len(nrow(scenarios))) {
    scenario <- scenarios[index, ]
    for (seed in seed_splits[[split]]) {
      directory <- run_directory(scenario$scenario_id, seed)
      if (!file.exists(file.path(directory, "counts.tsv"))) {
        stop(
          "Missing simulation ", directory,
          ". Run this script with --simulate first.",
          call. = FALSE
        )
      }
      run_metrics <- tryCatch(fit_one_run(directory), error = function(e) {
        # A small library can leave fewer negative-control guides than
        # `bb_calibrate_controls()` requires. That limit predates the
        # moderation and stops BARCS-original just as it stops the moderated
        # variants, so the run is recorded as skipped rather than silently
        # dropped or partially reported.
        message("skipped ", basename(directory), ": ", conditionMessage(e))
        NULL
      })
      if (is.null(run_metrics)) {
        skipped <- c(skipped, basename(directory))
        next
      }
      run_metrics$split <- split
      run_metrics$scenario_id <- scenario$scenario_id
      run_metrics$seed <- seed
      run_metrics$replicates <- scenario$replicates
      run_metrics$moi <- scenario$moi
      run_metrics$high_quality_guide_fraction <-
        scenario$high_quality_guide_fraction
      run_metrics$genes_simulated <- scenario$genes
      metric_rows[[length(metric_rows) + 1L]] <- run_metrics
      message("evaluated ", split, " ", basename(directory))
    }
  }
}
metrics <- do.call(rbind, metric_rows)
if (length(skipped)) {
  message(
    "Skipped ", length(skipped), " run(s) for want of negative controls: ",
    paste(skipped, collapse = ", ")
  )
}
write.csv(
  metrics,
  file.path(
    "data", "derived", "crispulator_facs_improved_barcs_holdout_metrics.csv"
  ),
  row.names = FALSE
)

summary_metrics <- c(
  "auroc", "average_precision", "directional_recall", "realized_fdp", "f1",
  "inactive_gene_p_below_0_05"
)
summary_table <- do.call(rbind, lapply(names(seed_splits), function(split) {
  split_data <- metrics[metrics$split == split, , drop = FALSE]
  do.call(rbind, lapply(unique(split_data$method), function(method) {
    method_data <- split_data[split_data$method == method, , drop = FALSE]
    do.call(rbind, lapply(summary_metrics, function(metric) {
      values <- method_data[[metric]]
      data.frame(
        analysis_protocol = analysis_protocol, split = split,
        method = method, metric = metric,
        runs = length(values), mean = mean(values), sd = sd(values),
        se = sd(values) / sqrt(length(values))
      )
    }))
  }))
}))
write.csv(
  summary_table,
  file.path(
    "data", "derived", "crispulator_facs_improved_barcs_holdout_summary.csv"
  ),
  row.names = FALSE
)

reference <- "BARCS-moderated"
paired_table <- do.call(rbind, lapply(names(seed_splits), function(split) {
 split_data <- metrics[metrics$split == split, , drop = FALSE]
 do.call(rbind, lapply(
  setdiff(unique(split_data$method), reference), function(method) {
    do.call(rbind, lapply(summary_metrics, function(metric) {
      merged <- merge(
        split_data[split_data$method == method,
                   c("scenario_id", "seed", metric)],
        split_data[split_data$method == reference,
                   c("scenario_id", "seed", metric)],
        by = c("scenario_id", "seed"), suffixes = c("_method", "_reference")
      )
      difference <- merged[[paste0(metric, "_method")]] -
        merged[[paste0(metric, "_reference")]]
      half_width <- qt(0.975, df = length(difference) - 1L) *
        sd(difference) / sqrt(length(difference))
      data.frame(
        analysis_protocol = analysis_protocol, split = split, method = method,
        reference = reference, metric = metric, runs = length(difference),
        mean_difference = mean(difference), sd_difference = sd(difference),
        lower_95 = mean(difference) - half_width,
        upper_95 = mean(difference) + half_width
      )
    }))
  }
 ))
}))
write.csv(
  paired_table,
  file.path(
    "data", "derived", "crispulator_facs_improved_barcs_holdout_paired.csv"
  ),
  row.names = FALSE
)

cat("\nMeans by split:\n")
print(summary_table, row.names = FALSE)
cat("\nPaired differences against", reference, ":\n")
print(paired_table, row.names = FALSE)
