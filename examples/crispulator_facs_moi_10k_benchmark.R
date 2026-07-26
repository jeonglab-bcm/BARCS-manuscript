#!/usr/bin/env Rscript

# Genome-scale CRISPulator FACS benchmark: two BARCS gene statistics against
# the general count models, at two multiplicities of infection.
#
# This replaces the earlier 400-gene sensitivity grid for the headline
# simulation. The library is 10,000 genes and 50,000 guides, which is the
# order of a real genome-wide screen rather than a pilot, and only two BARCS
# modes are carried: BARCS-original, the historical signed-z statistic, and
# BARCS-moderated, the same statistic after guide-dispersion moderation.
#
# MOI 0.20 and 0.30 are both reported. MOI 0.20 is the main-text setting: it
# is the more common pooled-screen operating point, and the one where single
# integration per cell is the better approximation. MOI 0.30 is the
# supplementary setting.
#
# Stage one needs Julia and the pinned environment in `julia/`:
#
#     Rscript examples/crispulator_facs_moi_10k_benchmark.R --simulate
#
# Stage two fits and evaluates, and needs edgeR, limma, and DESeq2:
#
#     Rscript examples/crispulator_facs_moi_10k_benchmark.R
#
# MAGeCK-MLE is not included. It is not installable from PyPI or from source
# in this environment, so it could not be refitted at these settings; the
# 400-gene comparisons elsewhere in the repository retain it.

options(stringsAsFactors = FALSE)
analysis_protocol <- "crispulator-moi-10k-v1"
nominal_fdr <- 0.10
seeds <- c(20250724L, 20250725L, 20250726L)
moi_levels <- c(0.20, 0.30)
main_moi <- 0.20
n_genes <- 10000L
n_replicates <- 4L
guide_quality <- 0.90
simulation_root <- file.path("results", "crispulator_facs", "moi_10k")

compact_number <- function(value) {
  vapply(value, function(x) format(x, trim = TRUE, scientific = FALSE),
         character(1))
}
run_directory <- function(moi, seed) {
  file.path(simulation_root, sprintf(
    "moi%s_g%d_r%d__seed_%d", compact_number(moi), n_genes, n_replicates, seed
  ))
}

if ("--simulate" %in% commandArgs(trailingOnly = TRUE)) {
  for (moi in moi_levels) {
    for (seed in seeds) {
      target <- run_directory(moi, seed)
      if (file.exists(file.path(target, "counts.tsv"))) {
        next
      }
      status <- system2("julia", c(
        "--project=julia", file.path("julia", "simulate_crispulator_facs.jl"),
        target, n_replicates, seed, moi, guide_quality, n_genes
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
source(file.path("R", "bbreg.R"))

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

# BARCS on one design (a subset of sample types), returning both gene modes.
fit_barcs_design <- function(counts, sample_data, guide_truth, sample_types) {
  keep_sample <- sample_data$sample_type %in% sample_types
  design_data <- droplevels(sample_data[keep_sample, , drop = FALSE])
  y <- counts[, keep_sample, drop = FALSE]
  guide_result <- bb_screen(
    counts = y, totals = colSums(y), data = design_data,
    formula = ~ phenotype_z + replicate, term = "phenotype_z",
    guide = guide_truth$guide, gene = guide_truth$gene,
    min_total_count = 30,
    ncores = as.integer(Sys.getenv("BARCS_NCORES", "4"))
  )
  negative_control <- guide_truth$class == "negcontrol"
  original <- bb_calibrate_controls(guide_result, negative_control,
                                    alpha = 0.05)
  moderated_guides <- bb_moderate_dispersion(guide_result, trend = TRUE)
  moderated <- bb_calibrate_controls(
    moderated_guides, negative_control, alpha = 0.05, method = "qq_slope"
  )
  list(
    guide_result = guide_result,
    original = original,
    moderated = moderated,
    prior_df = attr(moderated_guides, "prior_df"),
    results = list(
      `BARCS-original` = bb_gene_original(original, min_guides = 1L),
      `BARCS-moderated` = bb_gene_original(moderated, min_guides = 1L)
    )
  )
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

  fit_start <- proc.time()
  low_bulk_high <- fit_barcs_design(
    counts, sample_data, guide_truth, c("low", "bulk", "high")
  )
  guide_seconds <- unname((proc.time() - fit_start)[["elapsed"]])
  # Two-tail ablation: does the unsorted bulk sample add anything?
  low_high <- fit_barcs_design(
    counts, sample_data, guide_truth, c("low", "high")
  )
  original <- low_bulk_high$original
  moderated <- low_bulk_high$moderated

  gene_results <- low_bulk_high$results
  names(low_high$results) <- paste0(names(low_high$results), " (low-high)")
  gene_results <- c(gene_results, low_high$results)

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
    negative_control_gene <- assessed$class == "negcontrol"
    data.frame(
      analysis_protocol = analysis_protocol,
      method = sub(" \\(low-high\\)$", "", method),
      design = if (grepl("low-high", method, fixed = TRUE)) {
        "low-high"
      } else {
        "low-bulk-high"
      },
      genes = nrow(assessed),
      active_genes = sum(active),
      auroc = auroc(active, score),
      average_precision = average_precision(active, score),
      effect_spearman_active = cor(
        assessed$estimate[active], assessed$theoretical_phenotype[active],
        method = "spearman"
      ),
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
      negative_control_gene_p_below_0_05 =
        mean(assessed$p_value[negative_control_gene] < 0.05),
      guide_fit_seconds = guide_seconds,
      prior_df = low_bulk_high$prior_df,
      control_scale_original = attr(original, "control_scale"),
      control_scale_moderated = attr(moderated, "control_scale")
    )
  }))
}

metric_rows <- list()
for (moi in moi_levels) {
  for (seed in seeds) {
    directory <- run_directory(moi, seed)
    if (!file.exists(file.path(directory, "counts.tsv"))) {
      stop(
        "Missing simulation ", directory,
        ". Run this script with --simulate first.",
        call. = FALSE
      )
    }
    run_metrics <- fit_one_run(directory)
    run_metrics$moi <- moi
    run_metrics$seed <- seed
    run_metrics$genes_simulated <- n_genes
    run_metrics$replicates <- n_replicates
    run_metrics$scope <- if (moi == main_moi) "main" else "supplementary"
    metric_rows[[length(metric_rows) + 1L]] <- run_metrics
    message("evaluated ", basename(directory))
  }
}
metrics <- do.call(rbind, metric_rows)
write.csv(
  metrics,
  file.path("data", "derived", "crispulator_facs_moi_10k_metrics.csv"),
  row.names = FALSE
)

summary_metrics <- c(
  "auroc", "average_precision", "effect_spearman_active", "calls",
  "directional_recall", "realized_fdp", "f1",
  "negative_control_gene_p_below_0_05"
)
summary_table <- do.call(rbind, lapply(moi_levels, function(moi) {
 do.call(rbind, lapply(unique(metrics$design), function(design) {
  moi_data <- metrics[metrics$moi == moi & metrics$design == design,
                      , drop = FALSE]
  do.call(rbind, lapply(unique(moi_data$method), function(method) {
    method_data <- moi_data[moi_data$method == method, , drop = FALSE]
    do.call(rbind, lapply(summary_metrics, function(metric) {
      values <- method_data[[metric]]
      data.frame(
        analysis_protocol = analysis_protocol, moi = moi, design = design,
        scope = method_data$scope[1L], method = method, metric = metric,
        runs = length(values), mean = mean(values), sd = sd(values),
        se = sd(values) / sqrt(length(values))
      )
    }))
  }))
 }))
}))
write.csv(
  summary_table,
  file.path("data", "derived", "crispulator_facs_moi_10k_summary.csv"),
  row.names = FALSE
)

reference <- "BARCS-moderated"
paired_table <- do.call(rbind, lapply(moi_levels, function(moi) {
 moi_data <- metrics[metrics$moi == moi & metrics$design == "low-bulk-high",
                     , drop = FALSE]
 do.call(rbind, lapply(
  setdiff(unique(moi_data$method), reference), function(method) {
    do.call(rbind, lapply(summary_metrics, function(metric) {
      merged <- merge(
        moi_data[moi_data$method == method, c("seed", metric)],
        moi_data[moi_data$method == reference, c("seed", metric)],
        by = "seed", suffixes = c("_method", "_reference")
      )
      difference <- merged[[paste0(metric, "_method")]] -
        merged[[paste0(metric, "_reference")]]
      half_width <- qt(0.975, df = length(difference) - 1L) *
        sd(difference) / sqrt(length(difference))
      data.frame(
        analysis_protocol = analysis_protocol, moi = moi, method = method,
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
  file.path("data", "derived", "crispulator_facs_moi_10k_paired.csv"),
  row.names = FALSE
)

cat("\nGenome-scale CRISPulator benchmark,", n_genes, "genes, gene FDR",
    nominal_fdr, "\n")
print(summary_table, row.names = FALSE)
cat("\nPaired differences against", reference, ":\n")
print(paired_table, row.names = FALSE)
