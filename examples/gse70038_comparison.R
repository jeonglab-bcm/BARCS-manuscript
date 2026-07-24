#!/usr/bin/env Rscript

# Head-to-head analysis of GEO GSE70038.
#
# The design follows Table 5 / Box 3 of the MAGeCKFlute protocol:
# every sample has an initial-condition intercept, Day 0 samples have no
# additional indicator, and each terminal sample has one binary cell-line
# condition. This script fits beta-binomial regression to all guides, combines
# directional guide-level evidence to genes for comparison, and reads results
# from the official MAGeCK-MLE executable.

source(file.path("R", "bbreg.R"))

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

jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}
top_depleted <- function(data, n = 200L) {
  head(data$gene[order(data$estimate, data$p_value)], n)
}

effect_pairs <- do.call(rbind, lapply(terms, function(term_name) {
  bb <- bb_gene[bb_gene$term == term_name, ]
  mg <- mageck_gene[mageck_gene$term == term_name, ]
  merged <- merge(
    bb[, c("gene", "estimate")],
    mg[, c("gene", "estimate")],
    by = "gene",
    suffixes = c("_bb", "_mageck")
  )
  merged$term <- term_name
  merged$bb_depletion_rank <- rank(
    merged$estimate_bb, ties.method = "average", na.last = "keep"
  )
  merged$mageck_depletion_rank <- rank(
    merged$estimate_mageck, ties.method = "average", na.last = "keep"
  )
  merged
}))
write.csv(
  effect_pairs,
  gzfile(file.path(result_dir, "effect_concordance_pairs.csv.gz")),
  row.names = FALSE
)

concordance <- do.call(rbind, lapply(terms, function(term_name) {
  merged <- effect_pairs[effect_pairs$term == term_name, ]
  bb <- bb_gene[bb_gene$term == term_name, ]
  mg <- mageck_gene[mageck_gene$term == term_name, ]
  data.frame(
    term = term_name,
    comparison = "beta-binomial vs MAGeCK-MLE",
    spearman = cor(
      merged$bb_depletion_rank, merged$mageck_depletion_rank,
      method = "pearson", use = "complete.obs"
    ),
    top200_depleted_jaccard = jaccard(
      top_depleted(bb), top_depleted(mg)
    )
  )
}))
write.csv(
  concordance,
  file.path(result_dir, "method_concordance.csv"),
  row.names = FALSE
)

validated_gene <- c("PKMYT1", "HEATR1", "TFAP2C", "FBXO42", "HDAC2", "WEE1")
all_gene_results <- rbind(bb_gene, mageck_gene)
known_result <- all_gene_results[
  all_gene_results$gene %in% validated_gene,
  c("model", "gene", "term", "n_guides", "estimate", "p_value", "fdr")
]
known_result$depletion_rank <- ave(
  all_gene_results$estimate,
  interaction(all_gene_results$model, all_gene_results$term),
  FUN = function(value) rank(value, ties.method = "min", na.last = "keep")
)[all_gene_results$gene %in% validated_gene]
write.csv(
  known_result,
  file.path(result_dir, "published_gene_results.csv"),
  row.names = FALSE
)

pdf(file.path("figures", "gse70038_head_to_head.pdf"),
    width = 10, height = 8, onefile = TRUE)
old_par <- par(mfrow = c(2, 2), mar = c(4.2, 4.4, 2.4, 1))
for (term_name in terms) {
  bb <- bb_gene[bb_gene$term == term_name, ]
  mg <- mageck_gene[mageck_gene$term == term_name, ]
  pair_bm <- merge(bb, mg, by = "gene", suffixes = c("_bb", "_mageck"))
  x_limit <- quantile(
    pair_bm$estimate_bb, c(0.005, 0.995), na.rm = TRUE
  )
  y_limit <- quantile(
    pair_bm$estimate_mageck, c(0.005, 0.995), na.rm = TRUE
  )
  plot(
    pair_bm$estimate_bb, pair_bm$estimate_mageck,
    pch = 16, cex = 0.35, col = "#00000025",
    xlab = "Beta-binomial median guide effect",
    ylab = "MAGeCK-MLE beta", main = term_name,
    xlim = x_limit, ylim = y_limit
  )
  abline(0, 1, lty = 2)
  special <- pair_bm$gene %in% validated_gene
  points(pair_bm$estimate_bb[special], pair_bm$estimate_mageck[special],
         pch = 16, col = "#D55E00")
  text(pair_bm$estimate_bb[special], pair_bm$estimate_mageck[special],
       pair_bm$gene[special], pos = 3, cex = 0.65)

  bb_rank <- rank(
    pair_bm$estimate_bb, ties.method = "min", na.last = "keep"
  )
  mg_rank <- rank(
    pair_bm$estimate_mageck, ties.method = "min", na.last = "keep"
  )
  plot(
    bb_rank, mg_rank,
    pch = 16, cex = 0.35, col = "#00000025",
    xlab = "Beta-binomial depletion rank",
    ylab = "MAGeCK-MLE depletion rank",
    main = "Rank concordance"
  )
  abline(0, 1, lty = 2)

  bb_sig <- -log10(pmax(bb$fdr, .Machine$double.xmin))
  mg_sig <- -log10(pmax(mg$fdr, .Machine$double.xmin))
  plot(
    stats::ecdf(bb_sig[is.finite(bb_sig)]),
    verticals = TRUE, do.points = FALSE,
    xlab = expression(-log[10](FDR)),
    ylab = "Empirical cumulative fraction",
    main = "Gene-level significance",
    col = "#0072B2", lwd = 2
  )
  lines(
    stats::ecdf(mg_sig[is.finite(mg_sig)]),
    verticals = TRUE, do.points = FALSE,
    col = "#D55E00", lwd = 2
  )
  legend(
    "bottomright",
    legend = c("beta-binomial", "MAGeCK-MLE"),
    col = c("#0072B2", "#D55E00"), lwd = 2, bty = "n"
  )

  known <- known_result[
    known_result$term == term_name &
      known_result$gene %in% validated_gene,
  ]
  effect_matrix <- xtabs(estimate ~ gene + model, data = known)
  barplot(
    t(effect_matrix),
    beside = TRUE,
    las = 2,
    ylab = "Estimated depletion effect",
    main = "Published validation genes",
    col = c("#0072B2", "#D55E00")
  )
  abline(h = 0, lty = 2)
  legend(
    "topright", legend = colnames(effect_matrix),
    fill = c("#0072B2", "#D55E00"),
    bty = "n", cex = 0.7
  )
}
par(old_par)
dev.off()

cat("\nTable-5-style design:\n")
print(design, row.names = FALSE)
cat("\nMethod concordance:\n")
print(concordance, row.names = FALSE, digits = 3)
cat("\nPublished validation genes:\n")
print(known_result, row.names = FALSE, digits = 3)
