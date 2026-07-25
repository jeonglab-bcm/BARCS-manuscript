#!/usr/bin/env Rscript

# Generate the committed R reference used by the dependency-free browser tests.
# Run from the BARCS repository root:
#
#   Rscript web/scripts/generate-r-parity-reference.R guides
#   Rscript web/scripts/generate-r-parity-reference.R genes
#
# The script writes CSV to stdout so fixture changes remain reviewable.

arguments <- commandArgs(trailingOnly = TRUE)
mode <- if (length(arguments)) arguments[[1L]] else "guides"
if (!mode %in% c("guides", "genes")) {
  stop("Mode must be `guides` or `genes`.", call. = FALSE)
}

source(file.path("R", "bbreg.R"))

count_table <- read.csv(
  file.path("web", "public", "example-counts.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
sample_data <- read.csv(
  file.path("web", "public", "example-metadata.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
sample_data$batch <- factor(sample_data$batch)
sample_names <- sample_data$sample
count_matrix <- as.matrix(count_table[, sample_names, drop = FALSE])
storage.mode(count_matrix) <- "double"
totals <- colSums(count_matrix)
control <- tolower(count_table$control) %in% c(
  "1", "true", "yes", "y", "control", "ntc"
)

specifications <- list(
  time = list(formula = ~ time, term = "time"),
  time_batch = list(formula = ~ time + batch, term = "time"),
  time_by_batch = list(
    formula = ~ time * batch,
    term = "time:batchB"
  ),
  batch = list(formula = ~ batch, term = "batchB")
)

fit_specification <- function(name, specification) {
  raw <- bb_screen(
    count_matrix,
    sample_data,
    specification$formula,
    specification$term,
    totals = totals,
    guide = count_table$guide,
    gene = count_table$gene
  )
  calibrated <- bb_calibrate_controls(raw, control)
  data.frame(
    model = name,
    term = specification$term,
    guide = raw$guide,
    gene = raw$gene,
    control = control,
    estimate = raw$estimate,
    raw_std_error = raw$std_error,
    raw_t_value = raw$t_value,
    df = raw$df,
    raw_p_value = raw$p_value,
    raw_fdr = raw$fdr,
    rho = raw$rho,
    converged = raw$converged,
    calibrated_std_error = calibrated$std_error,
    calibrated_t_value = calibrated$t_value,
    calibrated_p_value = calibrated$p_value,
    calibrated_fdr = calibrated$fdr,
    control_scale = attr(calibrated, "control_scale"),
    stringsAsFactors = FALSE
  )
}

guide_reference <- do.call(
  rbind,
  Map(fit_specification, names(specifications), specifications)
)
rownames(guide_reference) <- NULL

if (mode == "guides") {
  write.csv(guide_reference, stdout(), row.names = FALSE, quote = TRUE)
} else {
  selected <- guide_reference$model == "time_batch"
  guide_result <- data.frame(
    gene = guide_reference$gene[selected],
    guide = guide_reference$guide[selected],
    estimate = guide_reference$estimate[selected],
    std_error = guide_reference$calibrated_std_error[selected],
    raw_std_error = guide_reference$raw_std_error[selected],
    t_value = guide_reference$calibrated_t_value[selected],
    df = guide_reference$df[selected],
    p_value = guide_reference$calibrated_p_value[selected],
    fdr = guide_reference$calibrated_fdr[selected],
    converged = guide_reference$converged[selected],
    stringsAsFactors = FALSE
  )
  gene_result <- bb_gene_consistency(
    guide_result,
    control = control
  )
  write.csv(gene_result, stdout(), row.names = FALSE, quote = TRUE)
}
