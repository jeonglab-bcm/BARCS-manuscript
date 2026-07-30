#!/usr/bin/env Rscript

# Multivariable BARCS and MAGeCK-MLE analysis of GEO GSE70038.
#
# The design follows Table 5 / Box 3 of the MAGeCKFlute protocol:
# every sample has an initial-condition intercept, Day 0 samples have no
# additional indicator, and each terminal sample has one binary cell-line
# condition. This script fits beta-binomial regression to all guides, combines
# directional guide-level evidence to genes, and reads results from the
# official MAGeCK-MLE executable. Functional comparison of the complementary
# hit sets is performed by `gse70038_complementary_enrichment.R`.

source(file.path("R", "bbreg.R"))
source(file.path("R", "method_palette.R"))

geo_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/",
  "GSE70038/suppl/GSE70038_rawReadCounts.txt.gz"
)
raw_path <- file.path("data", "raw", "GSE70038_rawReadCounts.txt.gz")
derived_dir <- file.path("data", "derived")
result_dir <- file.path("results", "gse70038")
dir.create(dirname(raw_path), recursive = TRUE, showWarnings = FALSE)
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

if (!file.exists(raw_path)) {
  download.file(geo_url, raw_path, mode = "wb", quiet = FALSE)
}

raw <- read.delim(
  gzfile(raw_path),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
sample_names <- names(raw)[10:25]
counts <- as.matrix(raw[, sample_names])
storage.mode(counts) <- "double"
guide <- raw[["sgRNA.ID"]]
gene <- vapply(
  strsplit(sub("^/", "", raw[["Gene name"]]), "/", fixed = TRUE),
  function(value) value[1L],
  character(1L)
)
gene[is.na(gene) | gene == ""] <- "unannotated"
rownames(counts) <- guide

condition <- c(
  "Initial_condition", "Initial_condition", "GSC0131_end", "GSC0131_end",
  "Initial_condition", "Initial_condition", "GSC0827_end", "GSC0827_end",
  "Initial_condition", "Initial_condition", "NSCCB660_end", "NSCCB660_end",
  "Initial_condition", "Initial_condition", "NSCU5_end", "NSCU5_end"
)
terms <- c("GSC0131_end", "GSC0827_end", "NSCCB660_end", "NSCU5_end")
design <- data.frame(
  samples = sample_names,
  Initial_condition = 1L,
  GSC0131_end = as.integer(condition == "GSC0131_end"),
  GSC0827_end = as.integer(condition == "GSC0827_end"),
  NSCCB660_end = as.integer(condition == "NSCCB660_end"),
  NSCU5_end = as.integer(condition == "NSCU5_end"),
  check.names = FALSE
)
model_data <- design[, terms, drop = FALSE]
model_formula <- ~ GSC0131_end + GSC0827_end + NSCCB660_end + NSCU5_end

design_path <- file.path(derived_dir, "GSE70038_table5_design.tsv")
mageck_count_path <- file.path(derived_dir, "GSE70038_mageck_counts.tsv")
write.table(
  design, design_path,
  sep = "\t", quote = FALSE, row.names = FALSE
)
if (!file.exists(mageck_count_path)) {
  write.table(
    data.frame(sgRNA = guide, Gene = gene, counts, check.names = FALSE),
    mageck_count_path,
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

mageck_executable <- file.path(".venv", "bin", "mageck")
mageck_prefix <- file.path(result_dir, "mageck")
mageck_gene_path <- paste0(mageck_prefix, ".gene_summary.txt")
if (file.exists(mageck_executable) && !file.exists(mageck_gene_path)) {
  status <- system2(
    mageck_executable,
    c(
      "mle",
      "-k", mageck_count_path,
      "-d", design_path,
      "-n", mageck_prefix,
      "--norm-method", "median",
      "--permutation-round", "1",
      "--threads", "4"
    )
  )
  if (status != 0) {
    warning("MAGeCK-MLE exited with status ", status)
  }
}
if (!file.exists(mageck_gene_path)) {
  stop(
    "MAGeCK-MLE output is unavailable. Install the official executable at ",
    "`.venv/bin/mageck` and rerun this script.",
    call. = FALSE
  )
}

cores <- max(1L, min(4L, parallel::detectCores(logical = FALSE)))
guide_blocks <- split(
  seq_len(nrow(counts)),
  cut(seq_len(nrow(counts)), breaks = cores, labels = FALSE)
)

run_blocks <- function(block_function) {
  if (.Platform$OS.type == "unix" && cores > 1L) {
    pieces <- parallel::mclapply(guide_blocks, block_function, mc.cores = cores)
  } else {
    pieces <- lapply(guide_blocks, block_function)
  }
  do.call(rbind, pieces)
}

message("Fitting beta-binomial regression to ", nrow(counts), " guides ...")
library_total <- colSums(counts)
bb_matrix <- run_blocks(function(index) {
  output <- matrix(
    NA_real_,
    nrow = length(index),
    ncol = length(terms) * 3L + 2L
  )
  for (row in seq_along(index)) {
    g <- index[row]
    fit <- tryCatch(
      bbreg(counts[g, ], library_total, model_formula, model_data),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    output[row, seq_along(terms)] <-
      fit$coefficient_table[terms, "estimate"]
    output[row, length(terms) + seq_along(terms)] <-
      fit$coefficient_table[terms, "std_error"]
    output[row, 2L * length(terms) + seq_along(terms)] <-
      fit$coefficient_table[terms, "p_value"]
    output[row, 3L * length(terms) + 1L] <- fit$rho
    output[row, 3L * length(terms) + 2L] <- as.numeric(fit$converged)
  }
  rownames(output) <- index
  output
})
bb_matrix <- bb_matrix[order(as.integer(rownames(bb_matrix))), , drop = FALSE]

make_guide_long <- function(model, fit_matrix, extra_name, extra_value) {
  result <- do.call(rbind, lapply(seq_along(terms), function(j) {
    data.frame(
      model = model,
      guide = guide,
      gene = gene,
      term = terms[j],
      estimate = fit_matrix[, j],
      std_error = fit_matrix[, length(terms) + j],
      p_value = fit_matrix[, 2L * length(terms) + j],
      stringsAsFactors = FALSE
    )
  }))
  result$fdr <- ave(
    result$p_value,
    result$term,
    FUN = function(p) p.adjust(p, method = "BH")
  )
  result[[extra_name]] <- rep(extra_value, times = length(terms))
  result
}

bb_guide <- make_guide_long(
  "beta-binomial",
  bb_matrix,
  "dispersion",
  bb_matrix[, 3L * length(terms) + 1L]
)
write.csv(
  bb_guide,
  gzfile(file.path(result_dir, "beta_binomial_guide_results.csv.gz")),
  row.names = FALSE
)
gene_index <- split(seq_along(gene), gene)
aggregate_gene <- function(guide_result) {
  do.call(rbind, lapply(terms, function(term_name) {
    rows <- guide_result[guide_result$term == term_name, ]
    estimate <- vapply(gene_index, function(index) {
      median(rows$estimate[index], na.rm = TRUE)
    }, numeric(1L))
    p_value <- vapply(gene_index, function(index) {
      keep <- is.finite(rows$p_value[index]) & is.finite(rows$estimate[index])
      if (!any(keep)) return(NA_real_)
      signed_z <- sign(rows$estimate[index][keep]) *
        qnorm(pmax(rows$p_value[index][keep] / 2, .Machine$double.xmin),
              lower.tail = FALSE)
      combined_z <- sum(signed_z) / sqrt(length(signed_z))
      2 * pnorm(-abs(combined_z))
    }, numeric(1L))
    data.frame(
      model = guide_result$model[1L],
      gene = names(gene_index),
      term = term_name,
      n_guides = lengths(gene_index),
      estimate = estimate,
      p_value = p_value,
      fdr = p.adjust(p_value, method = "BH"),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }))
}

bb_gene <- aggregate_gene(bb_guide)
write.csv(bb_gene, file.path(result_dir, "beta_binomial_gene_results.csv"),
          row.names = FALSE)

mageck_gene <- NULL
if (file.exists(mageck_gene_path)) {
  mageck_raw <- read.delim(
    mageck_gene_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  mageck_gene <- do.call(rbind, lapply(terms, function(term_name) {
    beta_column <- paste0(term_name, "|beta")
    p_column <- paste0(term_name, "|wald-p-value")
    fdr_column <- paste0(term_name, "|wald-fdr")
    data.frame(
      model = "MAGeCK-MLE 0.5.9.5",
      gene = mageck_raw[[1L]],
      term = term_name,
      n_guides = mageck_raw[[2L]],
      estimate = mageck_raw[[beta_column]],
      p_value = mageck_raw[[p_column]],
      fdr = mageck_raw[[fdr_column]],
      stringsAsFactors = FALSE
    )
  }))
  write.csv(
    mageck_gene,
    file.path(result_dir, "mageck_gene_results.csv"),
    row.names = FALSE
  )
}

cat("\nTable-5-style design:\n")
print(design, row.names = FALSE)
cat(
  "\nSaved guide- and gene-level BARCS and MAGeCK-MLE results for ",
  length(terms), " terminal coefficients.\n",
  sep = ""
)
