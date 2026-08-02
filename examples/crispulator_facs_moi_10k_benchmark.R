#!/usr/bin/env Rscript

# Genome-scale CRISPulator FACS benchmark: two BARCS gene statistics against
# the CRISPR-specific gene callers, at two multiplicities of infection.
#
# The comparators are MAGeCK-MLE and CRISPhieRmix. Both were designed for
# pooled screens, which is the comparison a screen analyst actually faces.
# The general RNA-seq count models are deliberately not scored here: DESeq2
# still runs, but only to supply the guide-level log2 fold changes that
# CRISPhieRmix documents as its input, not as a competing gene caller.
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
# Stage two fits and evaluates, and needs DESeq2 and CRISPhieRmix:
#
#     Rscript examples/crispulator_facs_moi_10k_benchmark.R
#
# All methods are scored on the genes every method returned a finite result
# for, because MAGeCK-MLE drops genes its own filters reject and an unequal
# gene universe would make the metrics incomparable.
#
# MAGeCK-MLE uses the official 0.5.9.5 binary at `.venv/bin/mageck`, matching
# the rest of the repository. Its marker score is affinely mapped to zero-one
# because the MAGeCK initializer requires a reference design row whose
# non-intercept entries are all zero.

options(stringsAsFactors = FALSE)
analysis_protocol <- "crispulator-moi-10k-v1"
nominal_fdr <- 0.10
# Threshold scan, log spaced over five orders of magnitude. A genome-wide
# screen is not run at a gene FDR of 0.2, so the informative range is the
# strict end: with 10,000 genes the question is how each method behaves when
# asked for a handful of confident hits. The grid stops at 0.20, which is
# already past any operating point a screen would use.
thresholds <- c(
  1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 0.01, 0.03, 0.05, 0.10, 0.20
)
# Comparing F1 at a matched *nominal* cutoff compares methods sitting at
# different real error rates, so the ranking it produces depends on the cutoff.
# These are the realized-FDP levels at which F1 is compared instead.
matched_fdp_targets <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.10)
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

needed <- c("DESeq2", "CRISPhieRmix")
absent <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(absent)) {
  stop(
    "This comparison needs ", paste(absent, collapse = ", "),
    ". DESeq2 comes from Bioconductor, or on Debian/Ubuntu with ",
    "apt-get install r-bioc-deseq2. CRISPhieRmix is installed from source ",
    "with R CMD INSTALL after cloning github.com/timydaley/CRISPhieRmix; it ",
    "needs sn and nloptr, which in turn need gfortran and BLAS/LAPACK ",
    "headers (apt-get install gfortran libblas-dev liblapack-dev).",
    call. = FALSE
  )
}
suppressPackageStartupMessages({
  library(DESeq2)
  library(CRISPhieRmix)
})
source(file.path("R", "load_barcs.R"))

mageck <- file.path(".venv", "bin", "mageck")
if (!file.exists(mageck)) {
  stop("Official MAGeCK 0.5.9.5 is required at `.venv/bin/mageck`.",
       call. = FALSE)
}
mageck_environment <- paste0(
  "PATH=", normalizePath(dirname(mageck)), ":", Sys.getenv("PATH")
)

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

  # MAGeCK-MLE, official binary, on the same low-bulk-high samples.
  mageck_dir <- file.path(directory, "mageck")
  dir.create(mageck_dir, showWarnings = FALSE, recursive = TRUE)
  count_path <- file.path(mageck_dir, "counts.txt")
  design_path <- file.path(mageck_dir, "design.txt")
  control_path <- file.path(mageck_dir, "control_sgrna.txt")
  mle_prefix <- file.path(mageck_dir, "mle")
  mle_path <- paste0(mle_prefix, ".gene_summary.txt")
  if (!file.exists(mle_path)) {
    write.table(
      data.frame(
        sgRNA = guide_truth$guide[keep_guide],
        Gene = guide_truth$gene[keep_guide],
        round(y[keep_guide, , drop = FALSE]),
        check.names = FALSE
      ),
      count_path, sep = "\t", quote = FALSE, row.names = FALSE
    )
    writeLines(
      guide_truth$guide[keep_guide & guide_truth$class == "negcontrol"],
      control_path
    )
    # Affine map of the ordered marker onto zero-one, so the low samples of
    # the reference replicate form the all-zero reference row MAGeCK needs.
    marker <- design_data$phenotype_z
    marker <- (marker - min(marker)) / (max(marker) - min(marker))
    replicate_levels <- levels(design_data$replicate)
    mageck_design <- data.frame(
      samples = design_data$sample, baseline = 1, phenotype = marker
    )
    for (level in replicate_levels[-1L]) {
      mageck_design[[paste0("rep_", level)]] <-
        as.integer(design_data$replicate == level)
    }
    reference_row <- which(
      mageck_design$phenotype == 0 &
        rowSums(mageck_design[, -(1:2), drop = FALSE]) == 0
    )[1L]
    stopifnot(is.finite(reference_row))
    mageck_design <- mageck_design[
      c(reference_row, setdiff(seq_len(nrow(mageck_design)), reference_row)),
    ]
    write.table(
      mageck_design, design_path,
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    status <- system2(mageck, c(
      "mle", "-k", count_path, "-d", design_path, "-n", mle_prefix,
      "--norm-method", "none", "--control-sgrna", control_path,
      "--permutation-round", "1", "--no-permutation-by-group",
      "--threads", as.character(as.integer(Sys.getenv("BARCS_NCORES", "4")))
    ), env = mageck_environment, stdout = FALSE, stderr = FALSE)
    if (status != 0L || !file.exists(mle_path)) {
      stop("MAGeCK-MLE failed for ", directory, call. = FALSE)
    }
  }
  mle_raw <- read.delim(mle_path, check.names = FALSE)
  gene_results[["MAGeCK-MLE"]] <- data.frame(
    gene = mle_raw$Gene,
    estimate = mle_raw[["phenotype|beta"]],
    p_value = mle_raw[["phenotype|p-value"]],
    fdr = mle_raw[["phenotype|fdr"]]
  )

  # DESeq2 supplies guide-level log2 fold changes. It is CRISPhieRmix's
  # documented input, not a scored gene caller: nothing below reads its
  # p-values, only the shrunken effect per guide.
  deseq_result <- results(DESeq(DESeqDataSetFromMatrix(
    countData = round(y[keep_guide, , drop = FALSE]),
    colData = design_data, design = formula
  ), quiet = TRUE), name = "phenotype_z")

  # CRISPhieRmix: a hierarchical mixture over guides within a gene, with the
  # negative-control guides supplying the empirical null. BIMODAL = TRUE
  # because this screen has both phenotype-increasing and -decreasing genes;
  # the one-sided default would score only one of them. The control guides are
  # consumed as the null, so control genes get no gene-level call, exactly as
  # they get none from MAGeCK-MLE.
  guide_gene <- guide_truth$gene[keep_guide]
  guide_control <- guide_truth$class[keep_guide] == "negcontrol"
  log_fold_change <- deseq_result$log2FoldChange
  usable <- is.finite(log_fold_change)
  targeting <- usable & !guide_control
  set.seed(20250724L)
  mixture <- CRISPhieRmix::CRISPhieRmix(
    x = log_fold_change[targeting],
    geneIds = factor(guide_gene[targeting]),
    negCtrl = log_fold_change[usable & guide_control],
    BIMODAL = TRUE, VERBOSE = FALSE
  )
  mixture_genes <- as.character(mixture$genes)
  # CRISPhieRmix reports a local false-discovery rate rather than a p-value
  # and no signed effect, so the gene effect is the mean guide log fold change
  # and the local fdr stands in for the p-value wherever a ranking is needed.
  gene_effect <- tapply(
    log_fold_change[targeting], guide_gene[targeting], mean
  )
  gene_results[["CRISPhieRmix"]] <- data.frame(
    gene = mixture_genes,
    estimate = unname(gene_effect[mixture_genes]),
    p_value = mixture$locfdr,
    fdr = mixture$FDR
  )

  # MAGeCK-MLE drops genes its own filters reject, so without this every
  # method would be scored on a different gene universe and the metrics would
  # not be comparable. Restrict all methods to the genes every method
  # returned a finite result for.
  finite_genes <- lapply(gene_results, function(result) {
    usable <- is.finite(result$estimate) & is.finite(result$p_value) &
      is.finite(result$fdr)
    result$gene[usable]
  })
  common_genes <- Reduce(intersect, finite_genes)
  stopifnot(length(common_genes) > 0L)

  # Cache the per-gene results so a change of threshold grid costs nothing.
  gene_cache <- file.path(directory, "gene_results.rds")
  saveRDS(
    list(gene_results = gene_results, common_genes = common_genes,
         gene_truth = gene_truth),
    gene_cache
  )

  # Threshold scan on the same common gene set.
  threshold_rows <- do.call(rbind, lapply(names(gene_results), function(method) {
    scanned <- merge(
      gene_truth[gene_truth$gene %in% common_genes, , drop = FALSE],
      gene_results[[method]][, c("gene", "estimate", "p_value", "fdr")],
      by = "gene", all.x = TRUE
    )
    scanned <- scanned[
      is.finite(scanned$estimate) & is.finite(scanned$p_value) &
        is.finite(scanned$fdr), , drop = FALSE
    ]
    active <- scanned$active
    sign_match <- sign(scanned$estimate) == scanned$expected_sign
    do.call(rbind, lapply(thresholds, function(threshold) {
      called <- scanned$fdr < threshold
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
      data.frame(
        analysis_protocol = analysis_protocol,
        method = sub(" \\(low-high\\)$", "", method),
        design = if (grepl("low-high", method, fixed = TRUE)) {
          "low-high"
        } else {
          "low-bulk-high"
        },
        nominal_fdr = threshold,
        genes = nrow(scanned),
        active_genes = sum(active),
        calls = sum(called),
        precision = precision,
        recall = recall,
        directional_recall = mean(called[active] & sign_match[active]),
        realized_fdp = if (sum(called) > 0) {
          false_positive / sum(called)
        } else {
          0
        },
        f1 = if ((precision + recall) > 0) {
          2 * precision * recall / (precision + recall)
        } else {
          0
        }
      )
    }))
  }))
  attr(threshold_rows, "is_threshold_scan") <- TRUE

  point_rows <- do.call(rbind, lapply(names(gene_results), function(method) {
    # The negative-control diagnostic is a per-method property of the known
    # nulls, so it is taken from the method's own output rather than the
    # common set. MAGeCK-MLE omits control genes from its gene summary
    # entirely once their sgRNAs are declared, so it is NA by construction.
    own <- merge(
      gene_truth,
      gene_results[[method]][, c("gene", "p_value")],
      by = "gene"
    )
    own_control <- own$class == "negcontrol" & is.finite(own$p_value)
    negative_control_rate <- if (any(own_control)) {
      mean(own$p_value[own_control] < 0.05)
    } else {
      NA_real_
    }
    assessed <- merge(
      gene_truth[gene_truth$gene %in% common_genes, , drop = FALSE],
      gene_results[[method]][, c("gene", "estimate", "p_value", "fdr")],
      by = "gene", all.x = TRUE
    )
    assessed <- assessed[
      is.finite(assessed$estimate) & is.finite(assessed$p_value) &
        is.finite(assessed$fdr), , drop = FALSE
    ]
    stopifnot(nrow(assessed) == length(common_genes))
    score <- -log10(pmax(assessed$p_value, .Machine$double.xmin))
    called <- assessed$fdr < nominal_fdr
    active <- assessed$active
    sign_match <- sign(assessed$estimate) == assessed$expected_sign
    true_positive <- sum(called & active)
    false_positive <- sum(called & !active)
    false_negative <- sum(!called & active)
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
      negative_control_gene_p_below_0_05 = negative_control_rate,
      guide_fit_seconds = guide_seconds,
      prior_df = low_bulk_high$prior_df,
      control_scale_original = attr(original, "control_scale"),
      control_scale_moderated = attr(moderated, "control_scale")
    )
  }))
  list(point = point_rows, scan = threshold_rows)
}

metric_rows <- list()
scan_rows <- list()
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
    fitted <- fit_one_run(directory)
    run_metrics <- fitted$point
    scan_metrics <- fitted$scan
    scan_metrics$moi <- moi
    scan_metrics$seed <- seed
    scan_metrics$scope <- if (moi == main_moi) "main" else "supplementary"
    scan_rows[[length(scan_rows) + 1L]] <- scan_metrics
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

scan_metrics <- do.call(rbind, scan_rows)
write.csv(
  scan_metrics,
  file.path("data", "derived", "crispulator_facs_moi_10k_f1_by_fdr.csv"),
  row.names = FALSE
)
# F1 at matched realized FDP, interpolated within each run before averaging.
interpolate_at <- function(realized, f1, target) {
  ordering <- order(realized)
  realized <- realized[ordering]
  f1 <- f1[ordering]
  if (target < min(realized) || target > max(realized)) {
    return(NA_real_)
  }
  stats::approx(realized, f1, xout = target, ties = "ordered")$y
}
matched_rows <- do.call(rbind, lapply(
  split(scan_metrics,
        list(scan_metrics$moi, scan_metrics$design, scan_metrics$method,
             scan_metrics$seed),
        drop = TRUE),
  function(group) {
    do.call(rbind, lapply(matched_fdp_targets, function(target) {
      data.frame(
        analysis_protocol = analysis_protocol,
        moi = group$moi[1L], design = group$design[1L],
        method = group$method[1L], seed = group$seed[1L],
        matched_fdp = target,
        f1 = interpolate_at(group$realized_fdp, group$f1, target)
      )
    }))
  }
))
matched_summary <- do.call(rbind, lapply(
  split(matched_rows,
        list(matched_rows$moi, matched_rows$design, matched_rows$method,
             matched_rows$matched_fdp),
        drop = TRUE),
  function(group) {
    usable <- group$f1[is.finite(group$f1)]
    data.frame(
      analysis_protocol = analysis_protocol,
      moi = group$moi[1L], design = group$design[1L],
      method = group$method[1L], matched_fdp = group$matched_fdp[1L],
      runs = length(usable),
      f1 = if (length(usable)) mean(usable) else NA_real_,
      f1_sd = if (length(usable) > 1L) sd(usable) else NA_real_
    )
  }
))
matched_summary <- matched_summary[
  order(matched_summary$moi, matched_summary$design, matched_summary$method,
        matched_summary$matched_fdp),
]
write.csv(
  matched_summary,
  file.path(
    "data", "derived", "crispulator_facs_moi_10k_f1_at_matched_fdp.csv"
  ),
  row.names = FALSE
)

scan_summary <- do.call(rbind, lapply(
  split(
    scan_metrics,
    list(scan_metrics$moi, scan_metrics$design, scan_metrics$method,
         scan_metrics$nominal_fdr),
    drop = TRUE
  ),
  function(group) {
    data.frame(
      analysis_protocol = analysis_protocol,
      moi = group$moi[1L], design = group$design[1L],
      method = group$method[1L], nominal_fdr = group$nominal_fdr[1L],
      runs = nrow(group),
      precision = mean(group$precision), recall = mean(group$recall),
      directional_recall = mean(group$directional_recall),
      realized_fdp = mean(group$realized_fdp),
      realized_fdp_sd = sd(group$realized_fdp),
      f1 = mean(group$f1), f1_sd = sd(group$f1),
      calls = mean(group$calls)
    )
  }
))
scan_summary <- scan_summary[
  order(scan_summary$moi, scan_summary$design, scan_summary$method,
        scan_summary$nominal_fdr),
]
write.csv(
  scan_summary,
  file.path(
    "data", "derived", "crispulator_facs_moi_10k_f1_by_fdr_summary.csv"
  ),
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
      values <- values[is.finite(values)]
      data.frame(
        analysis_protocol = analysis_protocol, moi = moi, design = design,
        scope = method_data$scope[1L], method = method, metric = metric,
        runs = length(values),
        mean = if (length(values)) mean(values) else NA_real_,
        sd = if (length(values) > 1L) sd(values) else NA_real_,
        se = if (length(values) > 1L) {
          sd(values) / sqrt(length(values))
        } else {
          NA_real_
        }
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


# ---- F1 and calibration across nominal FDR thresholds ----------------------
method_colours <- c(
  `BARCS-original` = "#56B4E9",
  `BARCS-moderated` = "#0072B2",
  `MAGeCK-MLE` = "#D55E00",
  CRISPhieRmix = "#7A3E9D"
)
method_pch <- c(16, 17, 15, 8)
names(method_pch) <- names(method_colours)
curve_methods <- names(method_colours)
x_positions <- -log10(thresholds)

# Label the decades as powers of ten rather than as strings of zeros; the
# intermediate 3x thresholds get an unlabelled tick so the grid stays visible.
labelled_thresholds <- c(0.20, 0.10, 0.01, 1e-3, 1e-4, 1e-5, 1e-6)
x_labels <- as.expression(lapply(labelled_thresholds, function(value) {
  exponent <- round(log10(value))
  if (isTRUE(all.equal(value, 10^exponent)) && exponent <= -3) {
    bquote(10^.(exponent))
  } else {
    format(value, trim = TRUE, scientific = FALSE)
  }
}))

draw_curve <- function(moi, metric, title, ylim) {
  panel <- scan_summary[
    scan_summary$moi == moi & scan_summary$design == "low-bulk-high",
    ,
    drop = FALSE
  ]
  plot(
    NA, xlim = range(x_positions), ylim = ylim, xaxt = "n",
    xlab = "Nominal gene FDR",
    ylab = if (metric == "f1") "F1 score" else "Realized FDP",
    main = title, bty = "l"
  )
  axis(1, at = x_positions, labels = FALSE, tcl = -0.25)
  axis(1, at = -log10(labelled_thresholds), labels = x_labels)
  if (metric == "realized_fdp") {
    # Realized equals nominal: a method on this line reports what it delivers.
    lines(x_positions, thresholds, lty = 2, lwd = 1.5, col = "#333333")
  }
  for (method in curve_methods) {
    rows <- panel[panel$method == method, , drop = FALSE]
    rows <- rows[match(thresholds, rows$nominal_fdr), , drop = FALSE]
    lines(
      x_positions, rows[[metric]], type = "o", pch = method_pch[[method]],
      lwd = 2, col = method_colours[[method]]
    )
  }
}

pdf(
  file.path("figures", "crispulator_facs_moi_10k_f1_by_fdr.pdf"),
  width = 10.5, height = 8, useDingbats = FALSE
)
layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE),
       heights = c(1, 1, 0.22))
par(mar = c(4.5, 4.3, 2.8, 1))
draw_curve(main_moi, "f1", "A  F1: MOI 0.20", c(0, 0.95))
draw_curve(main_moi, "realized_fdp", "B  Realized FDP: MOI 0.20", c(0, 0.30))
draw_curve(0.30, "f1", "C  F1: MOI 0.30", c(0, 0.95))
draw_curve(0.30, "realized_fdp", "D  Realized FDP: MOI 0.30", c(0, 0.30))
par(mar = c(0, 0, 0, 0))
plot.new()
legend(
  "center", legend = curve_methods,
  col = unname(method_colours[curve_methods]),
  pch = unname(method_pch[curve_methods]),
  lty = 1, lwd = 2, horiz = TRUE, bty = "n", cex = 0.78
)
dev.off()

cat("\nGenome-scale CRISPulator benchmark,", n_genes, "genes, gene FDR",
    nominal_fdr, "\n")
print(summary_table, row.names = FALSE)
cat("\nPaired differences against", reference, ":\n")
print(paired_table, row.names = FALSE)
