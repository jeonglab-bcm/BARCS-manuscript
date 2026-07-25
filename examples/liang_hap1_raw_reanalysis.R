#!/usr/bin/env Rscript

# BARCS reanalysis of the four deposited Liang et al. HAP1 endpoint FASTQs.
#
# Prerequisite:
#   bash scripts/run_liang_hap1_real_case.sh
#
# That command verifies each ENA FASTQ against its published MD5, applies the
# documented anchor trimming and Bowtie guide-alignment rules, and writes the
# raw guide-count matrix consumed here.

options(stringsAsFactors = FALSE)

input_dir <- file.path("results", "liang_hap1_raw")
count_path <- file.path(input_dir, "liang-hap1-raw-counts.csv")
metadata_path <- file.path(input_dir, "liang-hap1-raw-metadata.csv")
guide_path <- file.path("data", "raw", "liang_cas13", "guide_library.tsv")
expression_path <- file.path(
  "data", "raw", "liang_cas13", "lncrna_expression.tsv"
)
required <- c(count_path, metadata_path, guide_path, expression_path)
if (any(!file.exists(required))) {
  stop(
    "Missing raw-count inputs. Run ",
    "`bash scripts/run_liang_hap1_real_case.sh` first."
  )
}

source(file.path("R", "bbreg.R"))
counts_input <- read.csv(count_path, check.names = FALSE)
metadata <- read.csv(metadata_path, check.names = FALSE)
guide_library <- read.delim(guide_path, check.names = FALSE)
expression <- read.delim(expression_path, check.names = FALSE)

sample_names <- metadata$sample
if (!all(sample_names %in% names(counts_input))) {
  stop("Raw count columns do not match the sample metadata.")
}
counts <- as.matrix(counts_input[, sample_names, drop = FALSE])
storage.mode(counts) <- "integer"
if (anyNA(counts) || any(counts < 0L)) {
  stop("Raw guide counts must be non-negative integers.")
}
observed_totals <- colSums(counts)
if (!identical(as.integer(metadata$total), as.integer(observed_totals))) {
  stop("Metadata totals do not equal the complete raw count-column sums.")
}

sample_data <- data.frame(
  replicate = factor(metadata$replicate, levels = c("R1", "R2")),
  day14 = as.integer(metadata$day14)
)
control <- tolower(counts_input$control) == "true"

message(
  "Fitting ", nrow(counts), " guides from four verified HAP1 FASTQs ..."
)
started <- proc.time()
guide_result <- bb_screen(
  counts = counts,
  totals = observed_totals,
  data = sample_data,
  formula = ~ replicate + day14,
  term = "day14",
  guide = counts_input$guide,
  gene = counts_input$gene,
  min_total_count = 30L,
  ncores = min(4L, parallel::detectCores(logical = FALSE))
)
elapsed <- unname((proc.time() - started)[["elapsed"]])
guide_result <- bb_calibrate_controls(
  guide_result,
  control = control,
  alpha = 0.05
)
gene_result <- bb_gene_consistency(
  guide_result,
  control = control,
  min_guides = 3L,
  alpha = 0.05,
  min_control_genes = 10L
)

write.csv(
  guide_result,
  gzfile(file.path(input_dir, "liang-hap1-raw-barcs-guides.csv.gz")),
  row.names = FALSE
)
write.csv(
  gene_result,
  file.path(input_dir, "liang-hap1-raw-barcs-genes.csv"),
  row.names = FALSE
)

guide_annotation <- guide_library[
  match(counts_input$guide, guide_library$sgrna),
]
essential_genes <- unique(
  guide_annotation$gene[
    guide_annotation$target_group == "essential protein-coding gene"
  ]
)
hap1_expression <- expression[expression$cell_line == "HAP1", ]
unexpressed_lncRNA <- hap1_expression$gene[
  is.finite(hap1_expression$tpm) & hap1_expression$tpm == 0
]
tested_gene <- gene_result[!gene_result$control_gene, ]
discoveries <- tested_gene[
  is.finite(tested_gene$fdr) &
    tested_gene$fdr < 0.10 &
    tested_gene$estimate < 0,
]
essential_tested <- intersect(essential_genes, tested_gene$gene)
null_tested <- intersect(unexpressed_lncRNA, tested_gene$gene)

summary <- data.frame(
  metric = c(
    paste0("mapped_reads_", sample_names),
    "guides_in_library",
    "guides_with_finite_fit",
    "residual_degrees_of_freedom",
    "control_calibration_scale",
    "tested_noncontrol_genes",
    "depletion_discoveries_fdr_0.10",
    "essential_genes_tested",
    "essential_genes_recovered_fdr_0.10",
    "unexpressed_lncrna_tested",
    "unexpressed_lncrna_called_fdr_0.10",
    "elapsed_seconds"
  ),
  value = c(
    as.numeric(observed_totals),
    nrow(counts),
    sum(is.finite(guide_result$p_value)),
    unique(guide_result$df)[1L],
    attr(guide_result, "control_scale"),
    nrow(tested_gene),
    nrow(discoveries),
    length(essential_tested),
    sum(
      essential_tested %in% discoveries$gene
    ),
    length(null_tested),
    sum(null_tested %in% discoveries$gene),
    elapsed
  )
)
write.csv(
  summary,
  file.path(input_dir, "liang-hap1-raw-analysis-summary.csv"),
  row.names = FALSE
)
write.csv(
  summary[summary$metric != "elapsed_seconds", ],
  file.path(
    "data", "derived", "liang_hap1_raw_analysis_summary.csv"
  ),
  row.names = FALSE
)

message(
  "Raw HAP1 analysis complete: ", nrow(discoveries),
  " depleted genes at FDR < 0.10; results written to ", input_dir, "."
)
