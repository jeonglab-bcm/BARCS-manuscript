#!/usr/bin/env Rscript

# Head-to-head analysis of the Liang et al. transcriptome-scale Cas13 fitness
# screens (Cell Genomics, 2026; PRJNA1344834).
#
# Comparators:
#   1. BARCS beta-binomial regression;
#   2. official MAGeCK-MLE;
#   3. edgeR-QL, DESeq2, and limma-voom.
#
# All five methods estimate the same continuous time slope across days 0, 7,
# and 14. Endpoint RRA summaries are intentionally excluded because they do not
# estimate that longitudinal coefficient. The deposited count values are
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

raw_dir <- file.path("data", "raw", "liang_cas13")
result_dir <- file.path("results", "liang_cas13")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(
  "guide_library.tsv", "lncrna_expression.tsv"
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
cell_lines <- c("HAP1", "HEK293FT", "MDA-MB-231", "THP1")

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
    "^Day(0|7|14) Replicate [12] +\\(Count\\)$",
    names(processed),
    value = TRUE
  )
  if (length(count_columns) < 5L || length(count_columns) > 6L) {
    stop("Unexpected longitudinal count layout for ", cell_line)
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
    "^Day(0|7|14) Replicate.*$", "\\1", count_columns
  ))
  replicate <- as.integer(sub(
    "^Day(?:0|7|14) Replicate ([12]).*$", "\\1", count_columns,
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

source(file.path("R", "load_barcs.R"))

combine_guides <- function(result, retain_controls = TRUE) {
  gene_label <- result$gene
  ungrouped_control <- gene_label == "non-targeting"
  gene_label[ungrouped_control] <- paste0(
    "NT_single_", result$guide[ungrouped_control]
  )
  is_control <- startsWith(gene_label, "NT_")
  keep_gene <- retain_controls | !is_control
  indices <- split(which(keep_gene), gene_label[keep_gene])
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
      p_value = 2 * pnorm(
        -abs(sum(signed_z) / sqrt(length(signed_z)))
      )
    )
  }))
  combined$fdr <- p.adjust(combined$p_value, method = "BH")
  combined
}

assign_control_pseudogenes <- function(cell_guide, cell_line) {
  is_control <- cell_guide$target_group == "non-targeting"
  control_index <- which(is_control)[order(cell_guide$sgrna[is_control])]
  target_sizes <- as.integer(table(cell_guide$gene[!is_control]))
  target_sizes <- target_sizes[target_sizes > 0L]
  if (!length(control_index) || !length(target_sizes)) {
    stop("Could not construct aggregation-matched control pseudo-genes.")
  }

  set.seed(20260731L + match(cell_line, cell_lines))
  sampled_sizes <- sample(
    target_sizes,
    size = length(control_index),
    replace = TRUE
  )
  assignment <- rep(NA_character_, nrow(cell_guide))
  records <- list()
  position <- 1L
  group_index <- 0L
  while (position <= length(control_index)) {
    group_index <- group_index + 1L
    group_size <- min(
      sampled_sizes[group_index],
      length(control_index) - position + 1L
    )
    members <- control_index[position:(position + group_size - 1L)]
    group_name <- sprintf("NT_%s_%04d", gsub("[^A-Za-z0-9]", "", cell_line),
                          group_index)
    assignment[members] <- group_name
    records[[group_index]] <- data.frame(
      cell_line = cell_line,
      pseudo_gene = group_name,
      n_guides = group_size
    )
    position <- position + group_size
  }
  list(
    assignment = assignment,
    records = do.call(rbind, records)
  )
}

calibrate_gene_controls <- function(result, alpha = 0.05,
                                    min_controls = 20L,
                                    min_scale = 1) {
  is_control <- startsWith(result$gene, "NT_")
  signed_z <- sign(result$effect) * qnorm(
    pmax(result$p_value / 2, .Machine$double.xmin),
    lower.tail = FALSE
  )
  valid_control <- is_control & is.finite(signed_z)
  if (sum(valid_control) < min_controls) {
    stop("Too few finite non-targeting pseudo-genes for control scaling.")
  }
  empirical_cutoff <- unname(quantile(
    abs(signed_z[valid_control]),
    probs = 1 - alpha,
    type = 8
  ))
  scale <- max(
    min_scale,
    empirical_cutoff / qnorm(1 - alpha / 2)
  )
  result$raw_p_value <- result$p_value
  result$p_value <- 2 * pnorm(-abs(signed_z / scale))
  result$fdr <- p.adjust(result$p_value, method = "BH")
  result$control_scale <- scale
  result
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

has_complete_trajectories <- function(sample_data) {
  trajectory_table <- table(sample_data$replicate, sample_data$time_14)
  nrow(trajectory_table) >= 2L &&
    ncol(trajectory_table) >= 3L &&
    all(trajectory_table > 0L)
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

  formula <- if (has_complete_trajectories(sample_data)) {
    ~ replicate + time_14
  } else {
    ~ time_14
  }
  design <- model.matrix(formula, data = sample_data)
  coefficient <- match("time_14", colnames(design))
  if (is.na(coefficient)) {
    stop("The longitudinal time coefficient is absent from the design.")
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
    result <- result[
      match(cell_guide$sgrna, result$guide),
      ,
      drop = FALSE
    ]
    if (anyNA(result$guide)) {
      stop("Cached guide results do not match the Liang guide library.")
    }
    result$gene <- cell_guide$gene
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
      "^time_14($|_)", DESeq2::resultsNames(dds), value = TRUE
    )
    if (length(coefficient_name) != 1L) {
      stop("Could not identify one DESeq2 time coefficient.")
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
  control_groups <- assign_control_pseudogenes(cell_guide, cell_line)
  analysis_guide <- cell_guide
  is_control <- analysis_guide$target_group == "non-targeting"
  analysis_guide$gene[is_control] <- control_groups$assignment[is_control]
  audit <- input$audit
  audit <- audit[order(
    audit$day, audit$processed_replicate
  ), ]
  counts <- input$counts[, match(audit$sample, colnames(input$counts)), drop = FALSE]
  sample_data <- data.frame(
    time_14 = audit$day / 14,
    replicate = factor(audit$processed_replicate)
  )

  # Include the estimand in every cache name. This prevents an earlier
  # two-endpoint fit from being mistaken for the current three-time-point fit.
  cell_stem <- paste0(gsub("-", "_", cell_line), "_longitudinal")
  barcs_guide_path <- file.path(
    result_dir, paste0(cell_stem, "_barcs_guide.csv.gz")
  )
  if (!file.exists(barcs_guide_path) ||
      identical(Sys.getenv("RERUN_LIANG"), "1")) {
    formula <- if (has_complete_trajectories(sample_data)) {
      ~ replicate + time_14
    } else {
      ~ time_14
    }
    barcs_guide <- bb_screen(
      counts = counts,
      totals = colSums(counts),
      data = sample_data,
      formula = formula,
      term = "time_14",
      guide = analysis_guide$sgrna,
      gene = analysis_guide$gene,
      min_total_count = 30,
      ncores = min(4L, parallel::detectCores(logical = FALSE))
    )
    write.csv(barcs_guide, gzfile(barcs_guide_path), row.names = FALSE)
  } else {
    barcs_guide <- read.csv(gzfile(barcs_guide_path))
  }
  barcs_guide <- barcs_guide[
    match(analysis_guide$sgrna, barcs_guide$guide),
    ,
    drop = FALSE
  ]
  barcs_guide$gene <- analysis_guide$gene
  barcs_gene <- combine_guides(barcs_guide)
  write.csv(
    barcs_gene,
    file.path(result_dir, paste0(cell_stem, "_barcs_gene.csv")),
    row.names = FALSE
  )
  general_count <- fit_general_count_methods(
    counts, sample_data, analysis_guide, cell_stem
  )

  # Use the same aggregation-matched non-targeting pseudo-genes for every
  # method. This calibrates the native gene statistic, including the
  # guide-to-gene aggregation step, rather than transferring a one-guide scale
  # to multi-guide targets.
  mageck_gene <- analysis_guide$gene
  control_path <- file.path(
    result_dir, paste0(cell_stem, "_non_targeting_guides.txt")
  )
  writeLines(analysis_guide$sgrna[is_control], control_path)
  count_path <- file.path(result_dir, paste0(cell_stem, "_counts.tsv"))
  write.table(
    data.frame(
      sgRNA = analysis_guide$sgrna,
      Gene = mageck_gene,
      counts,
      check.names = FALSE
    ),
    count_path,
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  design_path <- file.path(result_dir, paste0(cell_stem, "_design.tsv"))
  if (has_complete_trajectories(sample_data)) {
    design <- data.frame(
      samples = audit$sample,
      baseline = 1,
      replicate2 = as.integer(audit$processed_replicate == 2L),
      time_14 = audit$day / 14
    )
  } else {
    design <- data.frame(
      samples = audit$sample,
      baseline = 1,
      time_14 = audit$day / 14
    )
  }
  write.table(
    design, design_path,
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  mle_prefix <- file.path(
    result_dir, paste0(cell_stem, "_aggregation_matched_mageck_mle")
  )
  mle_path <- paste0(mle_prefix, ".gene_summary.txt")
  if (!file.exists(mle_path) || identical(Sys.getenv("RERUN_LIANG"), "1")) {
    status <- system2(
      mageck,
      c(
        "mle", "-k", count_path, "-d", design_path,
        "-n", mle_prefix,
        "--norm-method", "none",
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
    effect = mle_raw[["time_14|beta"]],
    p_value = mle_raw[["time_14|wald-p-value"]],
    fdr = mle_raw[["time_14|wald-fdr"]]
  )

  add_method <- function(x, method) {
    x <- calibrate_gene_controls(x)
    x <- x[!startsWith(x$gene, "NT_"), , drop = FALSE]
    x$fdr <- p.adjust(x$p_value, method = "BH")
    if (!"n_guides" %in% names(x)) {
      x$n_guides <- NA_integer_
    }
    x$method <- method
    x$cell_line <- cell_line
    x$input_scale <- "rounded_normalized_ComBat_pseudocount"
    x$depletion_score <- -sign(x$effect) *
      -log10(pmax(x$p_value, .Machine$double.xmin))
    x[, c(
      "gene", "n_guides", "effect", "p_value", "fdr",
      "method", "cell_line", "input_scale", "depletion_score",
      "control_scale"
    )]
  }
  longitudinal_scores <- do.call(rbind, list(
    add_method(barcs_gene, "BARCS"),
    add_method(mle_gene, "MAGeCK-MLE"),
    add_method(general_count[["edgeR-QL"]], "edgeR-QL"),
    add_method(general_count[["DESeq2"]], "DESeq2"),
    add_method(general_count[["limma-voom"]], "limma-voom")
  ))

  endpoint_keep <- audit$day %in% c(0L, 14L)
  endpoint_audit <- audit[endpoint_keep, , drop = FALSE]
  endpoint_counts <- counts[
    , match(endpoint_audit$sample, colnames(counts)), drop = FALSE
  ]
  endpoint_data <- data.frame(
    time_14 = endpoint_audit$day / 14,
    replicate = factor(endpoint_audit$processed_replicate)
  )
  endpoint_stem <- paste0(gsub("-", "_", cell_line), "_endpoints")
  endpoint_path <- file.path(
    result_dir, paste0(endpoint_stem, "_barcs_guide.csv.gz")
  )
  if (!file.exists(endpoint_path) ||
      identical(Sys.getenv("RERUN_LIANG"), "1")) {
    endpoint_guide <- bb_screen(
      counts = endpoint_counts,
      totals = colSums(endpoint_counts),
      data = endpoint_data,
      formula = ~ replicate + time_14,
      term = "time_14",
      guide = analysis_guide$sgrna,
      gene = analysis_guide$gene,
      min_total_count = 30,
      ncores = min(4L, parallel::detectCores(logical = FALSE))
    )
    write.csv(endpoint_guide, gzfile(endpoint_path), row.names = FALSE)
  } else {
    endpoint_guide <- read.csv(gzfile(endpoint_path))
  }
  endpoint_guide <- endpoint_guide[
    match(analysis_guide$sgrna, endpoint_guide$guide),
    ,
    drop = FALSE
  ]
  endpoint_guide$gene <- analysis_guide$gene
  endpoint_gene <- combine_guides(endpoint_guide)
  endpoint_scores <- add_method(endpoint_gene, "BARCS endpoints")

  list(
    longitudinal = longitudinal_scores,
    ablation = rbind(
      longitudinal_scores[longitudinal_scores$method == "BARCS", ],
      endpoint_scores
    ),
    control_groups = control_groups$records
  )
}

cell_results <- lapply(cell_lines, run_cell_line)
scores <- do.call(rbind, lapply(cell_results, `[[`, "longitudinal"))
rownames(scores) <- NULL
ablation_scores <- do.call(rbind, lapply(cell_results, `[[`, "ablation"))
rownames(ablation_scores) <- NULL
control_group_records <- do.call(
  rbind,
  lapply(cell_results, `[[`, "control_groups")
)
write.csv(
  control_group_records,
  file.path(
    "data", "derived", "liang_cas13_control_pseudogene_sizes.csv"
  ),
  row.names = FALSE
)

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
ablation_scores <- merge(
  ablation_scores, gene_group, by = "gene", all.x = TRUE
)
ablation_scores <- merge(
  ablation_scores,
  expression,
  by = c("cell_line", "gene"),
  all.x = TRUE
)
ablation_scores$truth <- ifelse(
  ablation_scores$target_group == "essential protein-coding gene", 1L,
  ifelse(
    ablation_scores$target_group == "long non-coding RNA" &
      is.finite(ablation_scores$tpm) & ablation_scores$tpm == 0,
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

ablation_groups <- split(
  ablation_scores[!is.na(ablation_scores$truth), ],
  interaction(
    ablation_scores$cell_line[!is.na(ablation_scores$truth)],
    ablation_scores$method[!is.na(ablation_scores$truth)],
    drop = TRUE
  )
)
ablation_metrics <- do.call(rbind, lapply(ablation_groups, function(x) {
  data.frame(
    cell_line = x$cell_line[1L],
    analysis = ifelse(
      x$method[1L] == "BARCS",
      "Days 0, 7, and 14",
      "Days 0 and 14"
    ),
    average_precision = average_precision(x$truth, x$depletion_score),
    recall_at_5pct_null_fpr = recall_at_null_fpr(
      x$truth, x$depletion_score
    ),
    essential_fdr_0_10_recall = mean(
      x$fdr[x$truth == 1L] < 0.10 &
        x$effect[x$truth == 1L] < 0,
      na.rm = TRUE
    ),
    null_p_lt_0_05 = mean(x$p_value[x$truth == 0L] < 0.05)
  )
}))
rownames(ablation_metrics) <- NULL
write.csv(
  ablation_metrics,
  file.path("data", "derived", "liang_cas13_timepoint_ablation.csv"),
  row.names = FALSE
)
ablation_summary <- aggregate(
  ablation_metrics[, c(
    "average_precision", "recall_at_5pct_null_fpr",
    "essential_fdr_0_10_recall", "null_p_lt_0_05"
  )],
  by = list(analysis = ablation_metrics$analysis),
  FUN = mean,
  na.rm = TRUE
)
write.csv(
  ablation_summary,
  file.path("data", "derived", "liang_cas13_timepoint_ablation_summary.csv"),
  row.names = FALSE
)

control_scales <- unique(scores[, c(
  "cell_line", "method", "control_scale"
)])
control_scales <- control_scales[order(
  control_scales$cell_line, control_scales$method
), ]
write.csv(
  control_scales,
  file.path("data", "derived", "liang_cas13_control_scales.csv"),
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

cat("\nLiang Cas13 head-to-head macro-average:\n")
print(summary_metrics, row.names = FALSE)
cat("\nBARCS time-point ablation macro-average:\n")
print(ablation_summary, row.names = FALSE)
