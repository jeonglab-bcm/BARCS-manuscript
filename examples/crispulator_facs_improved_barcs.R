#!/usr/bin/env Rscript

# Guide-dispersion moderation for BARCS, measured against the general count
# models on the committed baseline CRISPulator realization.
#
# The head-to-head comparison showed BARCS trading ranking and recall for
# calibration: its realized false-discovery proportion sat far below the
# nominal threshold while edgeR-QL, DESeq2, and limma-voom ran two to three
# times above it. The per-gene case study located the cost. BARCS estimates one
# beta-binomial dispersion per guide from that guide's own residuals, so a
# handful of guides whose residuals happen to look quiet report standard errors
# their data do not support. Those guides drive the negative-control tail, and
# `bb_calibrate_controls()` answers by inflating the standard error of every
# guide in the screen.
#
# `bb_moderate_dispersion()` corrects those guides individually instead,
# shrinking each guide's variance inflation toward a library-wide trend. Once
# the quiet guides are fixed the blanket penalty is no longer needed, so most
# guides are tested more sharply even though the noisy ones are tested more
# strictly. A strictly conservative one-way variant is reported alongside.
#
# This script needs no files from the ignored `results/` tree.

options(stringsAsFactors = FALSE)
analysis_protocol <- "crispulator-improved-barcs-v1"
nominal_fdr <- 0.10

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

data_dir <- file.path("data", "derived", "crispulator_facs")
count_table <- read.delim(file.path(data_dir, "counts.tsv"), check.names = FALSE)
sample_data <- read.delim(
  file.path(data_dir, "sample_design.tsv"), check.names = FALSE
)
guide_truth <- read.delim(
  file.path(data_dir, "guide_truth.tsv"), check.names = FALSE
)
gene_truth <- read.delim(
  file.path(data_dir, "gene_truth.tsv"), check.names = FALSE
)
gene_truth$active <- tolower(as.character(gene_truth$active)) == "true"
counts <- as.matrix(count_table[, sample_data$sample, drop = FALSE])
storage.mode(counts) <- "double"
stopifnot(identical(count_table$guide, guide_truth$guide))

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
model_matrix <- model.matrix(~ phenotype_z + replicate, data = design_data)
keep_guide <- rowSums(y) >= 30

guide_result <- bb_screen(
  counts = y, totals = colSums(y), data = design_data,
  formula = ~ phenotype_z + replicate, term = "phenotype_z",
  guide = guide_truth$guide, gene = guide_truth$gene,
  min_total_count = 30,
  ncores = as.integer(Sys.getenv("BARCS_NCORES", "4"))
)
negative_control <- guide_truth$class == "negcontrol"

# Current BARCS.
current <- bb_calibrate_controls(guide_result, negative_control, alpha = 0.05)

# Proposed BARCS: moderate guide dispersion first, then estimate the residual
# null scale from the whole control quantile-quantile band rather than one
# tail order statistic.
moderated <- bb_moderate_dispersion(
  guide_result, trend = TRUE, one_way = FALSE, borrow_df = TRUE
)
improved <- bb_calibrate_controls(
  moderated, negative_control, alpha = 0.05, method = "qq_slope"
)

# Ablations, so the contribution of each half is visible.
moderation_only <- bb_calibrate_controls(
  moderated, negative_control, alpha = 0.05
)
qq_slope_only <- bb_calibrate_controls(
  guide_result, negative_control, alpha = 0.05, method = "qq_slope"
)
# Strictly conservative variant, for reference: never lowers a guide variance
# and never claims the prior degrees of freedom.
conservative <- bb_calibrate_controls(
  bb_moderate_dispersion(
    guide_result, trend = TRUE, one_way = TRUE, borrow_df = FALSE
  ),
  negative_control, alpha = 0.05, method = "qq_slope"
)

gene_results <- list(
  `BARCS-original` = bb_gene_original(current, min_guides = 1L),
  `BARCS-improved` = bb_gene_original(improved, min_guides = 1L),
  `BARCS-moderation-only` = bb_gene_original(moderation_only, min_guides = 1L),
  `BARCS-qq-slope-only` = bb_gene_original(qq_slope_only, min_guides = 1L),
  `BARCS-moderation-one-way` = bb_gene_original(conservative, min_guides = 1L)
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
  colData = design_data,
  design = ~ phenotype_z + replicate
), quiet = TRUE), name = "phenotype_z")
gene_results[["DESeq2"]] <- bb_gene_original(data.frame(
  gene = guide_truth$gene[keep_guide],
  estimate = deseq_result$log2FoldChange,
  p_value = deseq_result$pvalue
), min_guides = 1L)

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
evaluate <- function(method) {
  assessed <- merge(
    gene_truth, gene_results[[method]][, c("gene", "estimate", "p_value", "fdr")],
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
      2 * true_positive / (2 * true_positive + false_positive + false_negative)
    } else {
      0
    },
    inactive_gene_p_below_0_05 = mean(assessed$p_value[!active] < 0.05)
  )
}
metrics <- do.call(rbind, lapply(names(gene_results), evaluate))

scales <- data.frame(
  analysis_protocol = analysis_protocol,
  estimated_prior_df = attr(moderated, "prior_df"),
  control_scale_current = attr(current, "control_scale"),
  control_scale_improved = attr(improved, "control_scale"),
  control_scale_moderation_only = attr(moderation_only, "control_scale"),
  control_scale_qq_slope_only = attr(qq_slope_only, "control_scale")
)

write.csv(
  metrics,
  file.path("data", "derived", "crispulator_facs_improved_barcs_baseline.csv"),
  row.names = FALSE
)
write.csv(
  scales,
  file.path("data", "derived", "crispulator_facs_improved_barcs_scales.csv"),
  row.names = FALSE
)

cat("\nNull-scale diagnostics:\n")
print(scales, row.names = FALSE)
cat("\nBaseline realization, gene FDR", nominal_fdr, ":\n")
print(metrics[order(-metrics$f1), ], row.names = FALSE)
