#!/usr/bin/env Rscript

# Head-to-head analysis of the Liang et al. transcriptome-scale Cas13 fitness
# screens (Cell Genomics, 2026; PRJNA1344834).
#
# Comparators:
#   1. Liang RRA: deposited RobustRankAggreg v1.2.1 results;
#   2. MAGeCK-RRA: official `mageck test` on deposited processed counts;
#   3. MAGeCK-MLE: official `mageck mle` with replicate blocking;
#   4. BARCS: beta-binomial regression with the same replicate/day design;
#   5. edgeR-QL, DESeq2, and limma-voom with the same design.
#
# Liang RRA and MAGeCK-RRA are different algorithms.  The former ranks guide
# fold changes using the CRAN RobustRankAggreg package; the latter is MAGeCK's
# native guide-ranking/RRA workflow.  The deposited count values are
# fractional after normalization, ComBat correction, and outlier processing.
# BARCS correctly requires integer counts, so the normalized values are
# explicitly rounded to the nearest pseudo-count and that same matrix is given
# to every newly fitted method.  This is therefore a same-input processed-count
# sensitivity analysis, not a likelihood-faithful raw-count analysis.
#
# Published calls are a comparator, not ground truth.  Performance is assessed
# using the study's known-essential protein-coding controls as positives and
# lncRNAs with TPM == 0 in that cell line as biological null controls.

options(stringsAsFactors = FALSE)
source(file.path("R", "method_palette.R"))

raw_dir <- file.path("data", "raw", "liang_cas13")
result_dir <- file.path("results", "liang_cas13")
figure_path <- file.path("figures", "liang_cas13_benchmark.pdf")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(figure_path), recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(
  "guide_library.tsv", "lncrna_expression.tsv",
  "published_liang_rra.tsv"
)
missing_inputs <- required_inputs[
  !file.exists(file.path(raw_dir, required_inputs))
]
if (length(missing_inputs)) {
  stop(
    "Missing Liang inputs: ", paste(missing_inputs, collapse = ", "),
    ". Run `Rscript scripts/prepare_liang_cas13.R` first."
  )
}

guide <- read.delim(
  file.path(raw_dir, "guide_library.tsv"),
  check.names = FALSE
)
expression <- read.delim(
  file.path(raw_dir, "lncrna_expression.tsv"),
  check.names = FALSE
)
published <- read.delim(
  file.path(raw_dir, "published_liang_rra.tsv"),
  check.names = FALSE
)
cell_lines <- c("HAP1", "HEK293FT", "K562", "MDA-MB-231", "THP1")

load_processed_counts <- function(cell_line) {
  processed_path <- file.path(
    raw_dir,
    paste0(
      "published_processed_counts_",
      gsub("-", "_", cell_line),
      ".tsv"
    )
  )
  processed <- read.delim(processed_path, check.names = FALSE)
  count_columns <- grep(
    "^Day(0|14) Replicate [12] +\\(Count\\)$",
    names(processed),
    value = TRUE
  )
  if (length(count_columns) < 3L || length(count_columns) > 4L) {
    stop("Unexpected endpoint-count layout for ", cell_line)
  }
  guide_index <- match(processed$sgrna, guide$sgrna)
  keep <- !is.na(guide_index)
  counts <- as.matrix(processed[keep, count_columns, drop = FALSE])
  storage.mode(counts) <- "double"
  if (any(!is.finite(counts)) || any(counts < 0)) {
    stop("Processed counts must be finite and non-negative for ", cell_line)
  }
  mean_absolute_rounding <- mean(abs(counts - round(counts)))
  maximum_absolute_rounding <- max(abs(counts - round(counts)))
  counts <- round(counts)
  day <- as.integer(sub(
    "^Day(0|14) Replicate.*$", "\\1", count_columns
  ))
  replicate <- as.integer(sub(
    "^Day(?:0|14) Replicate ([12]).*$", "\\1", count_columns,
    perl = TRUE
  ))
  sample <- sprintf(
    "%s_Day%02d_R%d", cell_line, day, replicate
  )
  colnames(counts) <- sample
  list(
    counts = counts,
    guide_index = guide_index[keep],
    audit = data.frame(
      cell_line = cell_line,
      day = day,
      processed_replicate = replicate,
      sample = sample,
      processed_column = count_columns,
      input_scale = "rounded_normalized_ComBat_pseudocount",
      mean_absolute_rounding = mean_absolute_rounding,
      maximum_absolute_rounding = maximum_absolute_rounding
    )
  )
}

processed_inputs <- setNames(
  lapply(cell_lines, load_processed_counts),
  cell_lines
)
sample_audit <- do.call(rbind, lapply(processed_inputs, `[[`, "audit"))
write.csv(
  sample_audit,
  file.path(result_dir, "processed_count_input_audit.csv"),
  row.names = FALSE
)
write.csv(
  sample_audit,
  file.path("data", "derived", "liang_cas13_input_audit.csv"),
  row.names = FALSE
)

# Use the CB2 RcppArmadillo WLS kernel if its compiled shared library exists.
cpp_library <- file.path("CB2", "src", paste0("CB2", .Platform$dynlib.ext))
if (file.exists(cpp_library) &&
    requireNamespace("Rcpp", quietly = TRUE) &&
    requireNamespace("RcppArmadillo", quietly = TRUE)) {
  suppressPackageStartupMessages(library(Rcpp))
  suppressPackageStartupMessages(library(RcppArmadillo))
  if (!"CB2" %in% names(getLoadedDLLs())) {
    dyn.load(cpp_library)
  }
  source(file.path("CB2", "R", "RcppExports.R"))
}
source(file.path("R", "bbreg.R"))

combine_guides <- function(result) {
  keep_gene <- result$gene != "non-targeting"
  indices <- split(which(keep_gene), result$gene[keep_gene])
  combined <- do.call(rbind, lapply(names(indices), function(gene_name) {
    index <- indices[[gene_name]]
    valid <- is.finite(result$estimate[index]) &
      is.finite(result$p_value[index])
    if ("converged" %in% names(result)) {
      valid <- valid & !is.na(result$converged[index]) &
        result$converged[index]
    }
    index <- index[valid]
    if (!length(index)) {
      return(NULL)
    }
    signed_z <- sign(result$estimate[index]) *
      qnorm(
        pmax(result$p_value[index] / 2, .Machine$double.xmin),
        lower.tail = FALSE
      )
    data.frame(
      gene = gene_name,
      n_guides = length(index),
      effect = median(result$estimate[index]),
      p_value = 2 * pnorm(-abs(sum(signed_z) / sqrt(length(signed_z))))
    )
  }))
  combined$fdr <- p.adjust(combined$p_value, method = "BH")
  combined
}

make_guide_result <- function(cell_guide, estimate, p_value,
                              statistic = NA_real_) {
  estimate <- as.numeric(estimate)
  p_value <- as.numeric(p_value)
  statistic <- rep_len(as.numeric(statistic), length(estimate))
  data.frame(
    gene = cell_guide$gene,
    guide = cell_guide$sgrna,
    estimate = estimate,
    statistic = statistic,
    p_value = p_value,
    converged = is.finite(estimate) & is.finite(p_value)
  )
}

fit_general_count_methods <- function(counts, sample_data, cell_guide,
                                      cell_stem) {
  required <- c("edgeR", "DESeq2", "limma")
  unavailable <- required[
    !vapply(required, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(unavailable)) {
    stop(
      "Liang general-count comparisons require: ",
      paste(unavailable, collapse = ", ")
    )
  }

  formula <- if (nlevels(sample_data$replicate) >= 2L &&
                 nrow(sample_data) >= 4L) {
    ~ replicate + day14
  } else {
    ~ day14
  }
  design <- model.matrix(formula, data = sample_data)
  coefficient <- match("day14", colnames(design))
  if (is.na(coefficient)) {
    stop("The day-14 coefficient is absent from the design.")
  }

  # The deposited matrix is already normalized and ComBat-corrected. We add
  # only library-size offsets, not a second composition-normalization step.
  library_size <- colSums(counts)
  relative_library_size <- library_size /
    exp(mean(log(library_size)))

  fit_or_read <- function(method, fit) {
    path <- file.path(
      result_dir,
      paste0(cell_stem, "_", method, "_guide.csv.gz")
    )
    if (!file.exists(path) || identical(Sys.getenv("RERUN_LIANG"), "1")) {
      result <- fit()
      write.csv(result, gzfile(path), row.names = FALSE)
    } else {
      result <- read.csv(gzfile(path))
    }
    result
  }

  edger <- fit_or_read("edger_ql", function() {
    dge <- edgeR::DGEList(
      counts = counts,
      lib.size = library_size
    )
    dge$samples$norm.factors <- 1
    dge <- edgeR::estimateDisp(dge, design, robust = TRUE)
    fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
    test <- edgeR::glmQLFTest(fit, coef = coefficient)
    make_guide_result(
      cell_guide,
      estimate = test$table$logFC,
      p_value = test$table$PValue,
      statistic = sign(test$table$logFC) *
        sqrt(pmax(test$table$F, 0))
    )
  })

  deseq2 <- fit_or_read("deseq2", function() {
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = round(counts),
      colData = sample_data,
      design = formula
    )
    DESeq2::sizeFactors(dds) <- relative_library_size
    dds <- suppressMessages(DESeq2::DESeq(
      dds,
      test = "Wald",
      quiet = TRUE,
      minReplicatesForReplace = Inf
    ))
    coefficient_name <- grep(
      "^day14($|_)", DESeq2::resultsNames(dds), value = TRUE
    )
    if (length(coefficient_name) != 1L) {
      stop("Could not identify one DESeq2 day-14 coefficient.")
    }
    table <- DESeq2::results(
      dds,
      name = coefficient_name,
      independentFiltering = FALSE,
      cooksCutoff = FALSE
    )
    make_guide_result(
      cell_guide,
      estimate = table$log2FoldChange,
      p_value = table$pvalue,
      statistic = table$stat
    )
  })

  limma <- fit_or_read("limma_voom", function() {
    dge <- edgeR::DGEList(
      counts = counts,
      lib.size = library_size
    )
    dge$samples$norm.factors <- 1
    transformed <- limma::voom(dge, design, plot = FALSE)
    fit <- limma::eBayes(
      limma::lmFit(transformed, design),
      robust = TRUE
    )
    table <- limma::topTable(
      fit,
      coef = coefficient,
      number = Inf,
      sort.by = "none"
    )
    make_guide_result(
      cell_guide,
      estimate = table$logFC,
      p_value = table$P.Value,
      statistic = table$t
    )
  })

  list(
    `edgeR-QL` = combine_guides(edger),
    DESeq2 = combine_guides(deseq2),
    `limma-voom` = combine_guides(limma)
  )
}

mageck <- file.path(".venv", "bin", "mageck")
if (!file.exists(mageck)) {
  stop("Official MAGeCK 0.5.9.5 is required at `.venv/bin/mageck`.")
}
mageck_environment <- paste0(
  "PATH=", normalizePath(dirname(mageck)), ":", Sys.getenv("PATH")
)

run_cell_line <- function(cell_line) {
  message("Analysing ", cell_line, " ...")
  input <- processed_inputs[[cell_line]]
  cell_guide <- guide[input$guide_index, ]
  audit <- input$audit
  audit <- audit[order(
    audit$day, audit$processed_replicate
  ), ]
  counts <- input$counts[, match(audit$sample, colnames(input$counts)), drop = FALSE]
  sample_data <- data.frame(
    day14 = as.integer(audit$day == 14L),
    replicate = factor(audit$processed_replicate)
  )

  cell_stem <- gsub("-", "_", cell_line)
  barcs_guide_path <- file.path(
    result_dir, paste0(cell_stem, "_barcs_guide.csv.gz")
  )
  if (!file.exists(barcs_guide_path) ||
      identical(Sys.getenv("RERUN_LIANG"), "1")) {
    formula <- if (nlevels(sample_data$replicate) >= 2L &&
                   nrow(sample_data) >= 4L) {
      ~ replicate + day14
    } else {
      ~ day14
    }
    barcs_guide <- bb_screen(
      counts = counts,
      totals = colSums(counts),
      data = sample_data,
      formula = formula,
      term = "day14",
      guide = cell_guide$sgrna,
      gene = cell_guide$gene,
      min_total_count = 30,
      ncores = min(4L, parallel::detectCores(logical = FALSE))
    )
    write.csv(barcs_guide, gzfile(barcs_guide_path), row.names = FALSE)
  } else {
    barcs_guide <- read.csv(gzfile(barcs_guide_path))
  }
  barcs_calibrated_guide <- bb_calibrate_controls(
    barcs_guide,
    cell_guide$target_group == "non-targeting"
  )
  barcs_gene <- combine_guides(barcs_calibrated_guide)
  write.csv(
    barcs_gene,
    file.path(result_dir, paste0(cell_stem, "_barcs_gene.csv")),
    row.names = FALSE
  )
  general_count <- fit_general_count_methods(
    counts, sample_data, cell_guide, cell_stem
  )

  # Give each non-targeting guide its own pseudo-gene.  This prevents MAGeCK
  # from treating 1,000 unrelated controls as one enormous biological target.
  mageck_gene <- cell_guide$gene
  is_control <- cell_guide$target_group == "non-targeting"
  mageck_gene[is_control] <- paste0("NT_", cell_guide$sgrna[is_control])
  control_path <- file.path(
    result_dir, paste0(cell_stem, "_non_targeting_guides.txt")
  )
  writeLines(cell_guide$sgrna[is_control], control_path)
  count_path <- file.path(result_dir, paste0(cell_stem, "_counts.tsv"))
  write.table(
    data.frame(
      sgRNA = cell_guide$sgrna,
      Gene = mageck_gene,
      counts,
      check.names = FALSE
    ),
    count_path,
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  rra_prefix <- file.path(result_dir, paste0(cell_stem, "_mageck_rra"))
  rra_path <- paste0(rra_prefix, ".gene_summary.txt")
  if (!file.exists(rra_path) || identical(Sys.getenv("RERUN_LIANG"), "1")) {
    treatment <- paste(audit$sample[audit$day == 14L], collapse = ",")
    control <- paste(audit$sample[audit$day == 0L], collapse = ",")
    arguments <- c(
      "test", "-k", count_path,
      "-t", treatment, "-c", control,
      "--norm-method", "none",
      "--control-sgrna", control_path,
      "--gene-lfc-method", "median",
      "-n", rra_prefix
    )
    if (sum(audit$day == 0L) == sum(audit$day == 14L)) {
      arguments <- c(arguments, "--paired")
    }
    status <- system2(mageck, arguments, env = mageck_environment)
    if (status != 0 || !file.exists(rra_path)) {
      stop("MAGeCK-RRA failed for ", cell_line)
    }
  }
  rra_raw <- read.delim(rra_path, check.names = FALSE)
  rra_gene <- data.frame(
    gene = rra_raw$id,
    n_guides = rra_raw$num,
    effect = rra_raw[["neg|lfc"]],
    p_value = rra_raw[["neg|p-value"]],
    fdr = rra_raw[["neg|fdr"]]
  )

  design_path <- file.path(result_dir, paste0(cell_stem, "_design.tsv"))
  if (nlevels(sample_data$replicate) >= 2L && nrow(sample_data) >= 4L) {
    design <- data.frame(
      samples = audit$sample,
      baseline = 1,
      replicate2 = as.integer(audit$processed_replicate == 2L),
      day14 = as.integer(audit$day == 14L)
    )
  } else {
    design <- data.frame(
      samples = audit$sample,
      baseline = 1,
      day14 = as.integer(audit$day == 14L)
    )
  }
  write.table(
    design, design_path,
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  mle_prefix <- file.path(result_dir, paste0(cell_stem, "_mageck_mle"))
  mle_path <- paste0(mle_prefix, ".gene_summary.txt")
  if (!file.exists(mle_path) || identical(Sys.getenv("RERUN_LIANG"), "1")) {
    status <- system2(
      mageck,
      c(
        "mle", "-k", count_path, "-d", design_path,
        "-n", mle_prefix,
        "--norm-method", "none",
        "--control-sgrna", control_path,
        "--permutation-round", "1",
        "--no-permutation-by-group",
        "--threads", "4"
      ),
      env = mageck_environment
    )
    if (status != 0 || !file.exists(mle_path)) {
      stop("MAGeCK-MLE failed for ", cell_line)
    }
  }
  mle_raw <- read.delim(mle_path, check.names = FALSE)
  mle_gene <- data.frame(
    gene = mle_raw$Gene,
    n_guides = mle_raw$sgRNA,
    effect = mle_raw[["day14|beta"]],
    p_value = mle_raw[["day14|wald-p-value"]],
    fdr = mle_raw[["day14|wald-fdr"]]
  )

  liang <- published[published$cell_line == cell_line, ]
  liang_gene <- data.frame(
    gene = liang$gene,
    effect = liang$day14_log2_fold_change,
    p_value = liang$day14_p_value
  )
  liang_gene$fdr <- p.adjust(liang_gene$p_value, method = "BH")

  add_method <- function(x, method) {
    if (!"n_guides" %in% names(x)) {
      x$n_guides <- NA_integer_
    }
    x$method <- method
    x$cell_line <- cell_line
    x$input_scale <- "rounded_normalized_ComBat_pseudocount"
    x$depletion_score <- if (method %in% c(
      "Liang RRA", "MAGeCK-RRA"
    )) {
      -log10(pmax(x$p_value, .Machine$double.xmin))
    } else {
      -sign(x$effect) *
        -log10(pmax(x$p_value, .Machine$double.xmin))
    }
    x[, c(
      "gene", "n_guides", "effect", "p_value", "fdr",
      "method", "cell_line", "input_scale", "depletion_score"
    )]
  }
  do.call(rbind, list(
    add_method(liang_gene, "Liang RRA"),
    add_method(rra_gene, "MAGeCK-RRA"),
    add_method(mle_gene, "MAGeCK-MLE"),
    add_method(barcs_gene, "BARCS"),
    add_method(general_count[["edgeR-QL"]], "edgeR-QL"),
    add_method(general_count[["DESeq2"]], "DESeq2"),
    add_method(general_count[["limma-voom"]], "limma-voom")
  ))
}

scores <- do.call(rbind, lapply(cell_lines, run_cell_line))
rownames(scores) <- NULL

gene_group <- unique(guide[, c("gene", "target_group")])
group_priority <- match(
  gene_group$target_group,
  c(
    "essential protein-coding gene", "long non-coding RNA",
    "protein-coding gene", "non-targeting"
  )
)
gene_group <- gene_group[
  order(group_priority),
]
gene_group <- gene_group[!duplicated(gene_group$gene), ]
scores <- merge(scores, gene_group, by = "gene", all.x = TRUE)
scores <- merge(
  scores,
  expression,
  by = c("cell_line", "gene"),
  all.x = TRUE
)
scores$truth <- ifelse(
  scores$target_group == "essential protein-coding gene", 1L,
  ifelse(
    scores$target_group == "long non-coding RNA" &
      is.finite(scores$tpm) & scores$tpm == 0,
    0L,
    NA_integer_
  )
)
write.csv(
  scores,
  gzfile(file.path(result_dir, "all_gene_scores.csv.gz")),
  row.names = FALSE
)

roc_auc <- function(label, score) {
  keep <- is.finite(label) & is.finite(score)
  label <- label[keep]
  rank_value <- rank(score[keep], ties.method = "average")
  positives <- sum(label == 1L)
  negatives <- sum(label == 0L)
  (sum(rank_value[label == 1L]) -
     positives * (positives + 1) / 2) / (positives * negatives)
}

average_precision <- function(label, score) {
  keep <- is.finite(label) & is.finite(score)
  label <- label[keep][order(score[keep], decreasing = TRUE)]
  precision <- cumsum(label == 1L) / seq_along(label)
  mean(precision[label == 1L])
}

recall_at_null_fpr <- function(label, score, target_fpr = 0.05) {
  keep <- is.finite(label) & is.finite(score)
  label <- label[keep]
  score <- score[keep]
  cutoff <- unname(quantile(
    score[label == 0L],
    probs = 1 - target_fpr,
    type = 1,
    na.rm = TRUE
  ))
  mean(score[label == 1L] > cutoff)
}

metric_groups <- split(
  scores[!is.na(scores$truth), ],
  interaction(
    scores$cell_line[!is.na(scores$truth)],
    scores$method[!is.na(scores$truth)],
    drop = TRUE
  )
)
metrics <- do.call(rbind, lapply(metric_groups, function(x) {
  data.frame(
    cell_line = x$cell_line[1L],
    method = x$method[1L],
    positives = sum(x$truth == 1L),
    nulls = sum(x$truth == 0L),
    auroc = roc_auc(x$truth, x$depletion_score),
    average_precision = average_precision(x$truth, x$depletion_score),
    recall_at_5pct_null_fpr = recall_at_null_fpr(
      x$truth, x$depletion_score
    ),
    null_p_lt_0_05 = mean(x$p_value[x$truth == 0L] < 0.05),
    null_calibration_error = abs(
      mean(x$p_value[x$truth == 0L] < 0.05) - 0.05
    ),
    essential_fdr_0_10_recall = mean(
      x$fdr[x$truth == 1L] < 0.10 &
        x$effect[x$truth == 1L] < 0,
      na.rm = TRUE
    )
  )
}))
rownames(metrics) <- NULL
write.csv(
  metrics,
  file.path(result_dir, "benchmark_metrics_by_cell_line.csv"),
  row.names = FALSE
)
write.csv(
  metrics,
  file.path("data", "derived", "liang_cas13_metrics_by_cell_line.csv"),
  row.names = FALSE
)

summary_metrics <- aggregate(
  metrics[, c(
    "auroc", "average_precision", "recall_at_5pct_null_fpr",
    "null_p_lt_0_05", "null_calibration_error",
    "essential_fdr_0_10_recall"
  )],
  by = list(method = metrics$method),
  FUN = mean,
  na.rm = TRUE
)
write.csv(
  summary_metrics,
  file.path(result_dir, "benchmark_metrics_macro_average.csv"),
  row.names = FALSE
)
write.csv(
  summary_metrics,
  file.path("data", "derived", "liang_cas13_metrics_macro_average.csv"),
  row.names = FALSE
)

effect_wide <- reshape(
  scores[, c("cell_line", "gene", "method", "effect")],
  idvar = c("cell_line", "gene"),
  timevar = "method",
  direction = "wide"
)
effect_columns <- grep("^effect\\.", names(effect_wide), value = TRUE)
concordance <- do.call(rbind, lapply(cell_lines, function(cell_line) {
  x <- effect_wide[effect_wide$cell_line == cell_line, ]
  pairs <- combn(effect_columns, 2L, simplify = FALSE)
  do.call(rbind, lapply(pairs, function(pair) {
    data.frame(
      cell_line = cell_line,
      method_1 = sub("^effect\\.", "", pair[1L]),
      method_2 = sub("^effect\\.", "", pair[2L]),
      spearman = cor(
        x[[pair[1L]]], x[[pair[2L]]],
        method = "spearman", use = "complete.obs"
      )
    )
  }))
}))
write.csv(
  concordance,
  file.path(result_dir, "effect_rank_concordance.csv"),
  row.names = FALSE
)
write.csv(
  concordance,
  file.path("data", "derived", "liang_cas13_effect_rank_concordance.csv"),
  row.names = FALSE
)

method_order <- c(
  "Liang RRA", "MAGeCK-RRA", "MAGeCK-MLE", "BARCS",
  "edgeR-QL", "DESeq2", "limma-voom"
)
method_colour <- c(
  `Liang RRA` = barcs_method_colours[["Liang RRA"]],
  `MAGeCK-RRA` = barcs_method_colours[["MAGeCK-RRA"]],
  `MAGeCK-MLE` = barcs_method_colours[["MAGeCK"]],
  BARCS = barcs_method_colours[["BARCS"]],
  `edgeR-QL` = barcs_method_colours[["edgeR"]],
  DESeq2 = barcs_method_colours[["DESeq2"]],
  `limma-voom` = barcs_method_colours[["limma-voom"]]
)

pdf(figure_path, width = 11, height = 8.5)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
par(mar = c(7, 4, 3, 1))
for (metric in c(
  "average_precision", "recall_at_5pct_null_fpr",
  "null_calibration_error", "essential_fdr_0_10_recall"
)) {
  values <- summary_metrics[[metric]][
    match(method_order, summary_metrics$method)
  ]
  lower_better <- metric == "null_calibration_error"
  best <- if (lower_better) which.min(values) else which.max(values)
  bars <- barplot(
    values,
    names.arg = method_order,
    col = unname(method_colour[method_order]),
    border = ifelse(seq_along(values) == best, "black", NA),
    lwd = ifelse(seq_along(values) == best, 3, 1),
    las = 2,
    ylim = c(0, max(values, 0.05, na.rm = TRUE) * 1.15),
    ylab = switch(
      metric,
      average_precision = "Average precision",
      recall_at_5pct_null_fpr = "Essential recall at 5% null FPR",
      null_calibration_error = "|Null p < 0.05 fraction - 0.05|",
      essential_fdr_0_10_recall = "Essential recall at FDR 0.10"
    ),
    main = switch(
      metric,
      average_precision = "A. Essential-gene ranking",
      recall_at_5pct_null_fpr = "B. Equal null-FPR operating point",
      null_calibration_error = "C. Nominal null calibration",
      essential_fdr_0_10_recall = "D. Gene FDR threshold"
    )
  )
  text(
    bars, values,
    labels = sprintf("%.3f", values),
    pos = 3, cex = 0.8
  )
}
dev.off()

cat("\nLiang Cas13 head-to-head macro-average:\n")
print(summary_metrics, row.names = FALSE)
