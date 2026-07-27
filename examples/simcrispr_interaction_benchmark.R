#!/usr/bin/env Rscript

# Gene-by-environment interaction recovery on simCRISPR screens.
#
# CRISPulator generates a phenotype and sorts cells on it, which exercises the
# ordered-bin design but gives every guide a single effect. simCRISPR (Zhu et
# al., 2026) instead simulates a factorial screen: cells with and without
# knockout induction, treated and untreated, over several days, with PCR and
# sequencing noise applied on top. Each sgRNA therefore carries a knockout
# effect, a treatment effect, and an interaction between them.
#
# That interaction is the estimand BARCS was built for. It is a coefficient in
# a design matrix, not a group contrast, so it cannot be recovered by
# dichotomizing the experiment into two arms. This benchmark asks how well the
# interaction coefficient itself is recovered.
#
# Three properties of the simulator make the evaluation unusually clean:
#
#   * non-targeting guides have an interaction of exactly zero, so they are a
#     ground-truth null rather than an empirical approximation to one;
#   * safe-harbor guides also have zero targeted interaction but do carry a
#     cutting-related one, so they separate "no effect" from "no gene effect";
#   * the targeting guides have a continuous interaction, so effect recovery
#     can be scored without thresholding at all.
#
# The second axis is what BARCS uses as the beta-binomial denominator. Its
# default is the full-library total, on the argument that sequencing precision
# is set by depth. That argument fails when a large share of the library
# responds to the perturbation: dropout in the treated knockout arm inflates
# every surviving guide's share, and the whole screen acquires an apparent
# interaction. Normalizing the denominator on control guides instead removes
# it. Both control classes are tried, because they are not interchangeable:
# safe-harbor guides are cut and carry the DNA-damage response, non-targeting
# guides are not cut and do not. This is the comparison simCRISPR was built to
# make, and it is run here against an independent inference engine.
#
# Control guides used for normalization or calibration cannot also serve as the
# null, so each class is split in half. One half normalizes and calibrates and
# is declared to MAGeCK-MLE; the other half is held out and never seen, so the
# false-positive rate measured on it is not circular.
#
# Growth is logistic rather than exponential. Unbounded exponential growth over
# five days lets a single guide reach a fifth of the library and drives a
# 175-fold spread in library size, which is not what a real pooled screen looks
# like; the carrying capacity keeps both in a plausible range.
#
# CRISPhieRmix is not run here. It is defined by borrowing strength across the
# guides within a gene, and simCRISPR gives every sgRNA an independent
# interaction, so there is no gene structure for it to pool over. Grouping
# guides into synthetic genes would invent the very signal it exploits.
#
#     Rscript examples/simcrispr_interaction_benchmark.R
#
# MAGeCK-MLE uses the official 0.5.9.5 binary at `.venv/bin/mageck`, with each
# sgRNA as its own gene because the truth is guide-level.

options(stringsAsFactors = FALSE)
analysis_protocol <- "simcrispr-interaction-v1"
nominal_fdr <- 0.10
thresholds <- c(1e-4, 1e-3, 0.01, 0.05, 0.10, 0.20)
seeds <- c(20250724L, 20250725L, 20250726L)
n_guides <- 2000L
n_nontargeting <- 300L
n_safe_harbor <- 200L
replicates <- 3L
days <- 5L
reads_per_guide_sample <- 500L
# A guide counts as carrying a substantial interaction at this magnitude. The
# simulator gives every targeting guide a nonzero interaction, most of them far
# too small for any method to resolve, so a magnitude has to be named before
# thresholded metrics mean anything. Guides between zero and this value are
# scored for effect recovery but excluded from the detection metrics.
substantial_interaction <- 0.2
simulation_root <- file.path("results", "simcrispr", "interaction")

needed <- c("simCRISPR")
absent <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(absent)) {
  stop(
    "This benchmark needs ", paste(absent, collapse = ", "),
    ". Install it with R CMD INSTALL after cloning ",
    "github.com/bachergroup/simCRISPR.",
    call. = FALSE
  )
}
source(file.path("R", "bbreg.R"))

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

# ---- simulation ------------------------------------------------------------

simulate_one <- function(seed) {
  directory <- file.path(simulation_root, sprintf("seed_%d", seed))
  cache <- file.path(directory, "simulation.rds")
  if (file.exists(cache)) {
    return(readRDS(cache))
  }
  dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  set.seed(seed)
  raw <- simCRISPR::sim_crispr(
    method = "logit", samples = "independent",
    n_total = n_guides, n_ntgt = n_nontargeting, n_sfhb = n_safe_harbor,
    days = days, reps = replicates
  )
  # The full matrix carries the initial library plus all four arms; the
  # initial-library columns are not part of the factorial contrast.
  full <- raw$sim_full
  analysis_columns <- grep("^y0_", colnames(full), value = TRUE, invert = TRUE)
  molecules <- full[, analysis_columns, drop = FALSE]
  counts <- simCRISPR::seq_add(
    molecules,
    totalDepth = as.numeric(reads_per_guide_sample) *
      nrow(molecules) * ncol(molecules)
  )
  truth <- raw$true_eff_t
  truth$guide <- rownames(truth)
  truth$class <- sub("[0-9]+$", "", truth$guide)
  simulation <- list(counts = as.matrix(counts), truth = truth)
  saveRDS(simulation, cache)
  simulation
}

# The four arms are encoded as two binary factors, so the interaction is one
# model-matrix coefficient rather than a comparison of fitted group means.
sample_design <- function(sample_names) {
  data.frame(
    sample = sample_names,
    knockout = as.integer(grepl("^ko", sample_names)),
    treatment = as.integer(grepl("_trt_", sample_names))
  )
}

# The beta-binomial denominator. "library" is the full column total, BARCS's
# default. The control variants hold the chosen control class at a constant
# share of the denominator, which is what removes a composition shift, and are
# rescaled to the mean library size so the denominators stay interpretable as
# library sizes. They are rounded because a beta-binomial denominator counts
# reads.
effective_totals <- function(counts, control, mode) {
  library_total <- colSums(counts)
  if (identical(mode, "library")) {
    return(library_total)
  }
  control_total <- colSums(counts[control, , drop = FALSE])
  if (any(control_total <= 0)) {
    stop("A sample has no reads on the normalizing control guides.",
         call. = FALSE)
  }
  totals <- control_total / (mean(control_total) / mean(library_total))
  # The denominator must still dominate every count it is a denominator for.
  totals <- pmax(round(totals), apply(counts, 2L, max))
  totals
}

# ---- one realization -------------------------------------------------------

evaluate_one <- function(seed) {
  simulation <- simulate_one(seed)
  counts <- simulation$counts
  truth <- simulation$truth
  design_data <- sample_design(colnames(counts))
  stopifnot(
    length(unique(paste(design_data$knockout, design_data$treatment))) == 4L
  )

  truth <- truth[match(rownames(counts), truth$guide), , drop = FALSE]
  is_nontargeting <- truth$class == "ntgt"
  is_safe_harbor <- truth$class == "sfhb"
  is_targeting <- truth$class == "sg"

  # Each control class is split in half: one half normalizes and calibrates,
  # the other is held out as an untouched null. Both methods see the same
  # calibration half.
  set.seed(seed)
  halve <- function(flag) {
    index <- which(flag)
    chosen <- sample(index, floor(length(index) / 2))
    used <- logical(length(flag))
    used[chosen] <- TRUE
    used
  }
  calibrate_control <- halve(is_nontargeting)
  normalize_safe_harbor <- halve(is_safe_harbor)
  held_out_control <- is_nontargeting & !calibrate_control

  # BARCS twice over, under each choice of denominator. The calibration
  # controls are the same throughout, so the only thing that varies is what
  # the guide proportion is a proportion of.
  normalizations <- list(
    library = "library",
    `non-targeting` = calibrate_control,
    `safe-harbor` = normalize_safe_harbor
  )
  results <- list()
  guide_seconds <- NA_real_
  for (normalization in names(normalizations)) {
    control_set <- normalizations[[normalization]]
    totals <- effective_totals(counts, control_set, normalization)
    fit_start <- proc.time()
    guide_result <- bb_screen(
      counts = counts, totals = totals, data = design_data,
      formula = ~ knockout * treatment, term = "knockout:treatment",
      guide = rownames(counts), min_total_count = 30,
      ncores = as.integer(Sys.getenv("BARCS_NCORES", "4"))
    )
    if (identical(normalization, "library")) {
      guide_seconds <- unname((proc.time() - fit_start)[["elapsed"]])
    }
    control_flag <- calibrate_control[match(guide_result$guide, rownames(counts))]
    original <- bb_calibrate_controls(guide_result, control_flag, alpha = 0.05)
    moderated_guides <- bb_moderate_dispersion(guide_result, trend = TRUE)
    moderated <- bb_calibrate_controls(
      moderated_guides, control_flag, alpha = 0.05, method = "qq_slope"
    )
    # The empirical-null scale is the mechanism worth recording. A denominator
    # that imports a composition shift pushes every guide off zero, the
    # control calibration inflates the scale to absorb it, and the power that
    # inflation costs is what separates the normalizations below.
    results[[paste0("BARCS-original [", normalization, "]")]] <- data.frame(
      guide = original$guide, estimate = original$estimate,
      p_value = original$p_value, fdr = original$fdr,
      control_scale = attr(original, "control_scale")
    )
    results[[paste0("BARCS-moderated [", normalization, "]")]] <- data.frame(
      guide = moderated$guide, estimate = moderated$estimate,
      p_value = moderated$p_value, fdr = moderated$fdr,
      control_scale = attr(moderated, "control_scale")
    )
  }

  # MAGeCK-MLE with one sgRNA per gene, because the truth is guide-level.
  directory <- file.path(simulation_root, sprintf("seed_%d", seed))
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
        sgRNA = rownames(counts), Gene = rownames(counts),
        round(counts), check.names = FALSE
      ),
      count_path, sep = "\t", quote = FALSE, row.names = FALSE
    )
    writeLines(rownames(counts)[calibrate_control], control_path)
    mageck_design <- data.frame(
      samples = design_data$sample, baseline = 1L,
      knockout = design_data$knockout, treatment = design_data$treatment,
      interaction = design_data$knockout * design_data$treatment
    )
    # The MAGeCK initializer needs a reference row whose non-baseline entries
    # are all zero: the untreated, no-knockout arm.
    reference_row <- which(
      mageck_design$knockout == 0L & mageck_design$treatment == 0L
    )[1L]
    stopifnot(is.finite(reference_row))
    mageck_design <- mageck_design[
      c(reference_row, setdiff(seq_len(nrow(mageck_design)), reference_row)),
    ]
    write.table(
      mageck_design, design_path, sep = "\t", quote = FALSE, row.names = FALSE
    )
    status <- system2(mageck, c(
      "mle", "-k", count_path, "-d", design_path, "-n", mle_prefix,
      "--norm-method", "control", "--control-sgrna", control_path,
      "--permutation-round", "1", "--no-permutation-by-group",
      "--threads", as.character(as.integer(Sys.getenv("BARCS_NCORES", "4")))
    ), env = mageck_environment, stdout = FALSE, stderr = FALSE)
    if (status != 0L || !file.exists(mle_path)) {
      stop("MAGeCK-MLE failed for ", directory, call. = FALSE)
    }
  }
  mle_raw <- read.delim(mle_path, check.names = FALSE)
  results[["MAGeCK-MLE [non-targeting]"]] <- data.frame(
    guide = mle_raw$Gene,
    estimate = mle_raw[["interaction|beta"]],
    p_value = mle_raw[["interaction|p-value"]],
    fdr = mle_raw[["interaction|fdr"]],
    control_scale = NA_real_
  )

  # Score every method on the guides all of them returned a finite result for.
  finite_guides <- lapply(results, function(result) {
    result$guide[is.finite(result$estimate) & is.finite(result$p_value) &
                   is.finite(result$fdr)]
  })
  common_guides <- Reduce(intersect, finite_guides)
  stopifnot(length(common_guides) > 0L)

  truth_row <- match(common_guides, truth$guide)
  interaction_truth <- truth$INT[truth_row]
  class_common <- truth$class[truth_row]
  held_out_common <- held_out_control[match(common_guides, rownames(counts))]
  targeting_common <- class_common == "sg"
  safe_harbor_common <- class_common == "sfhb"

  # Detection labels: substantial targeting interactions are the positives,
  # held-out true zeros are the negatives, the rest is an ambiguous band.
  positive <- targeting_common & abs(interaction_truth) >= substantial_interaction
  negative <- held_out_common
  scored <- positive | negative

  point_rows <- list()
  scan_rows <- list()
  for (method in names(results)) {
    result <- results[[method]][match(common_guides, results[[method]]$guide), ]
    score <- -log10(pmax(result$p_value, .Machine$double.xmin))
    called <- result$fdr < nominal_fdr
    sign_match <- sign(result$estimate) == sign(interaction_truth)

    true_positive <- sum(called[scored] & positive[scored])
    false_positive <- sum(called[scored] & negative[scored])
    false_negative <- sum(!called[scored] & positive[scored])
    precision <- if ((true_positive + false_positive) > 0) {
      true_positive / (true_positive + false_positive)
    } else {
      0
    }
    recall <- true_positive / sum(positive)

    point_rows[[method]] <- data.frame(
      analysis_protocol = analysis_protocol,
      method = method,
      seed = seed,
      guides = length(common_guides),
      substantial_guides = sum(positive),
      held_out_controls = sum(negative),
      effect_spearman = cor(
        result$estimate[targeting_common], interaction_truth[targeting_common],
        method = "spearman"
      ),
      auroc = auroc(positive[scored], score[scored]),
      average_precision = average_precision(positive[scored], score[scored]),
      calls = sum(called[scored]),
      directional_recall = mean(
        called[positive] & sign_match[positive]
      ),
      realized_fdp = if (sum(called[scored]) > 0) {
        false_positive / sum(called[scored])
      } else {
        0
      },
      f1 = if ((precision + recall) > 0) {
        2 * precision * recall / (precision + recall)
      } else {
        0
      },
      # The assumption-free numbers: these guides have an interaction of
      # exactly zero and neither method was allowed to see them.
      held_out_control_call_rate = mean(called[negative]),
      safe_harbor_call_rate = mean(called[safe_harbor_common]),
      control_scale = result$control_scale[1L],
      guide_fit_seconds = guide_seconds
    )

    scan_rows[[method]] <- do.call(rbind, lapply(thresholds, function(cut) {
      cut_called <- result$fdr < cut
      cut_true <- sum(cut_called[scored] & positive[scored])
      cut_false <- sum(cut_called[scored] & negative[scored])
      cut_precision <- if ((cut_true + cut_false) > 0) {
        cut_true / (cut_true + cut_false)
      } else {
        0
      }
      cut_recall <- cut_true / sum(positive)
      data.frame(
        analysis_protocol = analysis_protocol, method = method, seed = seed,
        nominal_fdr = cut, calls = sum(cut_called[scored]),
        realized_fdp = if (sum(cut_called[scored]) > 0) {
          cut_false / sum(cut_called[scored])
        } else {
          0
        },
        held_out_control_call_rate = mean(cut_called[negative]),
        safe_harbor_call_rate = mean(cut_called[safe_harbor_common]),
        f1 = if ((cut_precision + cut_recall) > 0) {
          2 * cut_precision * cut_recall / (cut_precision + cut_recall)
        } else {
          0
        }
      )
    }))
  }
  list(
    point = do.call(rbind, point_rows),
    scan = do.call(rbind, scan_rows)
  )
}

# ---- run -------------------------------------------------------------------

point_list <- list()
scan_list <- list()
for (seed in seeds) {
  evaluated <- evaluate_one(seed)
  point_list[[length(point_list) + 1L]] <- evaluated$point
  scan_list[[length(scan_list) + 1L]] <- evaluated$scan
  message("evaluated simCRISPR seed ", seed)
}
metrics <- do.call(rbind, point_list)
scan <- do.call(rbind, scan_list)

dir.create(file.path("data", "derived"), showWarnings = FALSE, recursive = TRUE)
write.csv(
  metrics, file.path("data", "derived", "simcrispr_interaction_metrics.csv"),
  row.names = FALSE
)
write.csv(
  scan, file.path("data", "derived", "simcrispr_interaction_scan.csv"),
  row.names = FALSE
)

summarise <- function(frame, keys, value_columns) {
  do.call(rbind, lapply(split(frame, frame[keys], drop = TRUE), function(part) {
    do.call(rbind, lapply(value_columns, function(column) {
      values <- part[[column]]
      values <- values[is.finite(values)]
      if (!length(values)) {
        return(NULL)
      }
      cbind(
        unique(part[keys]),
        data.frame(
          metric = column, runs = length(values), mean = mean(values),
          sd = stats::sd(values),
          se = stats::sd(values) / sqrt(length(values))
        )
      )
    }))
  }))
}

metric_columns <- c(
  "effect_spearman", "auroc", "average_precision", "calls",
  "directional_recall", "realized_fdp", "f1",
  "held_out_control_call_rate", "safe_harbor_call_rate", "control_scale"
)
summary_table <- summarise(metrics, "method", metric_columns)
scan_summary <- summarise(
  scan, c("method", "nominal_fdr"),
  c("calls", "realized_fdp", "held_out_control_call_rate",
    "safe_harbor_call_rate", "f1")
)
write.csv(
  summary_table,
  file.path("data", "derived", "simcrispr_interaction_summary.csv"),
  row.names = FALSE
)
write.csv(
  scan_summary,
  file.path("data", "derived", "simcrispr_interaction_scan_summary.csv"),
  row.names = FALSE
)

cat("\nsimCRISPR interaction recovery,", n_guides, "guides,",
    replicates, "replicates per arm, over", length(seeds), "seeds\n")
wide <- reshape(
  summary_table[, c("method", "metric", "mean")],
  idvar = "method", timevar = "metric", direction = "wide"
)
names(wide) <- sub("^mean\\.", "", names(wide))
print(wide, row.names = FALSE)
cat("\nCalibration across cutoffs (held-out true-null call rate):\n")
null_rate <- scan_summary[
  scan_summary$metric == "held_out_control_call_rate",
  c("method", "nominal_fdr", "mean")
]
print(
  reshape(null_rate, idvar = "method", timevar = "nominal_fdr",
          direction = "wide"),
  row.names = FALSE
)
