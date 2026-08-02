#!/usr/bin/env Rscript

# Which genes do the general count models call that BARCS does not, and why?
#
# This runs on the committed baseline CRISPulator realization in
# `data/derived/crispulator_facs/`, so unlike the head-to-head script it needs
# no files from the ignored `results/` tree. The BARCS fit is checked against
# the committed head-to-head row for the same run before anything is reported.
#
# edgeR-QL, DESeq2, and limma-voom are fitted guide by guide on the same
# design and then aggregated with the same signed-z rule as BARCS-original, so
# the only thing that differs between the four methods is the guide-level
# count model.

options(stringsAsFactors = FALSE)
analysis_protocol <- "crispulator-miss-cases-v1"
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
source(file.path("R", "load_barcs.R"))

data_dir <- file.path("data", "derived", "crispulator_facs")
count_table <- read.delim(file.path(data_dir, "counts.tsv"), check.names = FALSE)
sample_data <- read.delim(
  file.path(data_dir, "sample_design.tsv"),
  check.names = FALSE
)
guide_truth <- read.delim(
  file.path(data_dir, "guide_truth.tsv"),
  check.names = FALSE
)
gene_truth <- read.delim(
  file.path(data_dir, "gene_truth.tsv"),
  check.names = FALSE
)
gene_truth$active <- tolower(as.character(gene_truth$active)) == "true"
counts <- as.matrix(count_table[, sample_data$sample, drop = FALSE])
storage.mode(counts) <- "double"
stopifnot(identical(count_table$guide, guide_truth$guide))

# Ordered-phenotype score: the conditional normal mean of each sorted bin.
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
min_total_count <- 30
keep_guide <- rowSums(y) >= min_total_count
model_matrix <- model.matrix(~ phenotype_z + replicate, data = design_data)
stopifnot(colnames(model_matrix)[2L] == "phenotype_z")

# ---- BARCS -----------------------------------------------------------------
guide_result <- bb_screen(
  counts = y,
  totals = colSums(y),
  data = design_data,
  formula = ~ phenotype_z + replicate,
  term = "phenotype_z",
  guide = guide_truth$guide,
  gene = guide_truth$gene,
  min_total_count = min_total_count,
  ncores = as.integer(Sys.getenv("BARCS_NCORES", "4"))
)
negative_control <- guide_truth$class == "negcontrol"
calibrated <- bb_calibrate_controls(
  guide_result,
  negative_control,
  alpha = 0.05
)
control_scale <- attr(calibrated, "control_scale")

# ---- general count models, same design, same gene combiner -----------------
dge <- DGEList(counts = y[keep_guide, , drop = FALSE])
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge, model_matrix, robust = TRUE)
quasi_likelihood <- glmQLFTest(
  glmQLFit(dge, model_matrix, robust = TRUE),
  coef = "phenotype_z"
)
edger_guide <- data.frame(
  gene = guide_truth$gene[keep_guide],
  estimate = quasi_likelihood$table$logFC,
  p_value = quasi_likelihood$table$PValue
)

voom_fit <- eBayes(lmFit(voom(dge, model_matrix), model_matrix))
limma_guide <- data.frame(
  gene = guide_truth$gene[keep_guide],
  estimate = voom_fit$coefficients[, "phenotype_z"],
  p_value = voom_fit$p.value[, "phenotype_z"]
)

deseq_dataset <- DESeqDataSetFromMatrix(
  countData = round(y[keep_guide, , drop = FALSE]),
  colData = design_data,
  design = ~ phenotype_z + replicate
)
deseq_result <- results(
  DESeq(deseq_dataset, quiet = TRUE),
  name = "phenotype_z"
)
deseq_guide <- data.frame(
  gene = guide_truth$gene[keep_guide],
  estimate = deseq_result$log2FoldChange,
  p_value = deseq_result$pvalue
)

gene_results <- list(
  `BARCS-original` = bb_gene_original(calibrated, min_guides = 1L),
  `edgeR-QL` = bb_gene_original(edger_guide, min_guides = 1L),
  DESeq2 = bb_gene_original(deseq_guide, min_guides = 1L),
  `limma-voom` = bb_gene_original(limma_guide, min_guides = 1L)
)

# ---- provenance check against the committed head-to-head row ---------------
committed <- read.csv(file.path(
  "data", "derived",
  "crispulator_facs_external_head_to_head_metrics.csv"
))
committed <- committed[
  committed$method == "BARCS-original" &
    committed$seed == 20250724L &
    committed$scenario_id == "moi_0.25_quality_0.9_genes_400_replicates_4",
]
stopifnot(nrow(committed) == 1L)

called_of <- function(method) {
  result <- gene_results[[method]]
  setNames(result$fdr < nominal_fdr, result$gene)[gene_truth$gene]
}
fdr_of <- function(method) {
  result <- gene_results[[method]]
  setNames(result$fdr, result$gene)[gene_truth$gene]
}
active <- gene_truth$active
barcs_called <- called_of("BARCS-original")
stopifnot(sum(barcs_called) == committed$calls_fdr_0_10)
message(
  "BARCS reproduces the committed run: ",
  sum(barcs_called), " calls, control scale ",
  format(round(control_scale, 4))
)

count_methods <- c("edgeR-QL", "DESeq2", "limma-voom")
count_called <- Reduce(`&`, lapply(count_methods, called_of))
missed <- active & count_called & !barcs_called

# ---- what the extra calls are worth ----------------------------------------
marginal_rows <- lapply(count_methods, function(method) {
  called <- called_of(method)
  extra <- called & !barcs_called
  data.frame(
    analysis_protocol = analysis_protocol,
    method = method,
    calls = sum(called),
    extra_over_barcs = sum(extra),
    extra_active = sum(extra & active),
    extra_precision = mean(active[extra]),
    barcs_only_calls = sum(barcs_called & !called)
  )
})
marginal_table <- do.call(rbind, marginal_rows)

# Extra calls stratified by how close BARCS came to calling them.
barcs_fdr <- fdr_of("BARCS-original")
edger_extra <- called_of("edgeR-QL") & !barcs_called
bands <- list(c(0.10, 0.15), c(0.15, 0.25), c(0.25, 0.50), c(0.50, 1.01))
band_table <- do.call(rbind, lapply(bands, function(band) {
  in_band <- edger_extra & barcs_fdr >= band[1L] & barcs_fdr < band[2L]
  data.frame(
    analysis_protocol = analysis_protocol,
    barcs_fdr_lower = band[1L],
    barcs_fdr_upper = band[2L],
    extra_calls = sum(in_band),
    extra_active = sum(in_band & active),
    extra_precision = if (sum(in_band) > 0) mean(active[in_band]) else NA_real_
  )
}))

# ---- per-gene table for the genes BARCS misses -----------------------------
calibrated$gene <- guide_truth$gene
mean_total <- mean(colSums(y))
variance_inflation <- 1 + (mean_total - 1) * calibrated$rho
finite_vif <- is.finite(variance_inflation)
gene_vif <- tapply(
  variance_inflation[finite_vif], calibrated$gene[finite_vif], mean
)
gene_knockdown <- tapply(guide_truth$knockdown, guide_truth$gene, mean)

miss_table <- data.frame(
  analysis_protocol = analysis_protocol,
  gene = gene_truth$gene[missed],
  class = gene_truth$class[missed],
  behavior = gene_truth$behavior[missed],
  theoretical_phenotype = gene_truth$theoretical_phenotype[missed],
  mean_knockdown = unname(gene_knockdown[gene_truth$gene[missed]]),
  mean_variance_inflation = unname(gene_vif[gene_truth$gene[missed]]),
  barcs_fdr = unname(barcs_fdr[missed]),
  edger_fdr = unname(fdr_of("edgeR-QL")[missed]),
  deseq_fdr = unname(fdr_of("DESeq2")[missed]),
  limma_fdr = unname(fdr_of("limma-voom")[missed])
)
miss_table <- miss_table[order(-abs(miss_table$theoretical_phenotype)), ]

# ---- how much of the gap is the one-way calibration ratchet? ---------------
sweep_table <- do.call(rbind, lapply(seq(1.00, 1.40, by = 0.04), function(s) {
  rescaled <- calibrated
  rescaled$p_value <- 2 * pt(
    -abs(calibrated$raw_t_value / s),
    df = calibrated$df
  )
  scaled_gene <- bb_gene_original(rescaled, min_guides = 1L)
  called <- setNames(
    scaled_gene$fdr < nominal_fdr, scaled_gene$gene
  )[gene_truth$gene]
  data.frame(
    analysis_protocol = analysis_protocol,
    control_scale = s,
    calls = sum(called),
    true_positives = sum(called & active),
    realized_fdp = mean(!active[called]),
    recovered_missed_genes = sum(called & missed)
  )
}))

write.csv(
  miss_table,
  file.path("data", "derived", "crispulator_facs_miss_cases.csv"),
  row.names = FALSE
)
write.csv(
  rbind(
    data.frame(marginal_table, row.names = NULL)
  ),
  file.path("data", "derived", "crispulator_facs_miss_marginal_value.csv"),
  row.names = FALSE
)
write.csv(
  band_table,
  file.path("data", "derived", "crispulator_facs_miss_precision_bands.csv"),
  row.names = FALSE
)
write.csv(
  sweep_table,
  file.path(
    "data", "derived", "crispulator_facs_calibration_scale_sweep.csv"
  ),
  row.names = FALSE
)

cat("\nMarginal value of the calls the count models add over BARCS:\n")
print(marginal_table, row.names = FALSE)
cat("\nExtra edgeR-QL calls by how close BARCS came:\n")
print(band_table, row.names = FALSE)
cat("\nActive genes all three count models call and BARCS misses:\n")
print(miss_table, row.names = FALSE)
cat("\nCalibration scale sweep (applied scale was ",
    format(round(control_scale, 4)), "):\n", sep = "")
print(sweep_table, row.names = FALSE)
