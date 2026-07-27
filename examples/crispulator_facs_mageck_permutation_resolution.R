#!/usr/bin/env Rscript

# Why MAGeCK-MLE cannot be compared at strict gene-FDR cutoffs.
#
# In the genome-scale scan, MAGeCK-MLE's call count is constant across several
# orders of magnitude of requested FDR. That is not a property of the data. Its
# gene p-value is a permutation tail probability estimated from `genes x rounds`
# draws, so it is quantized in steps of about 1/(genes x rounds), and every gene
# whose statistic exceeds all permuted values is reported with p exactly zero.
# Those genes are called at any cutoff, however strict, which produces a plateau.
#
# This script measures the effect directly by refitting one realization at
# several round counts. The prediction is that the quantization step and the
# number of genes at p = 0 both scale as 1/rounds, and that the plateau moves
# rather than disappears.
#
# It is slow: each additional round costs roughly one full gene-wise refit, so
# the 10-round fit takes about 80 minutes on four threads at 10,000 genes.
# Existing MAGeCK outputs are reused when present.
#
#     Rscript examples/crispulator_facs_mageck_permutation_resolution.R

options(stringsAsFactors = FALSE)
analysis_protocol <- "crispulator-mageck-permutation-resolution-v1"
rounds_grid <- c(1L, 10L)
cutoffs <- c(1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.05, 0.10)
run_directory <- file.path(
  "results", "crispulator_facs", "moi_10k",
  "moi0.2_g10000_r4__seed_20250724"
)
mageck_directory <- file.path(run_directory, "mageck")

source(file.path("R", "mageck.R"))
mageck <- mageck_executable()
mageck_check_version()

required <- file.path(
  mageck_directory, c("counts.txt", "design.txt", "control_sgrna.txt")
)
if (!all(file.exists(required))) {
  stop(
    "Run examples/crispulator_facs_moi_10k_benchmark.R first; it writes the ",
    "MAGeCK inputs this diagnostic reuses.",
    call. = FALSE
  )
}
mageck_environment <- paste0(
  "PATH=", normalizePath(dirname(mageck)), ":", Sys.getenv("PATH")
)

gene_truth <- read.delim(
  file.path(run_directory, "gene_truth.tsv"), check.names = FALSE
)
gene_truth$active <- tolower(as.character(gene_truth$active)) == "true"

summary_path <- function(rounds) {
  if (rounds == 1L) {
    file.path(mageck_directory, "mle.gene_summary.txt")
  } else {
    file.path(mageck_directory, sprintf("mle_p%d.gene_summary.txt", rounds))
  }
}

for (rounds in rounds_grid) {
  target <- summary_path(rounds)
  if (file.exists(target)) {
    next
  }
  prefix <- sub("\\.gene_summary\\.txt$", "", target)
  message("fitting MAGeCK with ", rounds, " permutation round(s); this is slow")
  status <- system2(mageck, c(
    "mle", "-k", file.path(mageck_directory, "counts.txt"),
    "-d", file.path(mageck_directory, "design.txt"), "-n", prefix,
    "--norm-method", "none",
    "--control-sgrna", file.path(mageck_directory, "control_sgrna.txt"),
    "--permutation-round", as.character(rounds),
    "--no-permutation-by-group",
    "--threads", as.character(as.integer(Sys.getenv("BARCS_NCORES", "4")))
  ), env = mageck_environment, stdout = FALSE, stderr = FALSE)
  if (status != 0L || !file.exists(target)) {
    stop("MAGeCK failed at ", rounds, " rounds.", call. = FALSE)
  }
}

rows <- list()
for (rounds in rounds_grid) {
  result <- read.delim(summary_path(rounds), check.names = FALSE)
  assessed <- merge(
    gene_truth, result[, c("Gene", "phenotype|p-value", "phenotype|fdr")],
    by.x = "gene", by.y = "Gene"
  )
  p_value <- assessed[["phenotype|p-value"]]
  fdr <- assessed[["phenotype|fdr"]]
  positive <- p_value[p_value > 0]
  for (cutoff in cutoffs) {
    called <- fdr < cutoff
    true_positive <- sum(called & assessed$active)
    false_positive <- sum(called & !assessed$active)
    precision <- if (sum(called) > 0) true_positive / sum(called) else 0
    recall <- true_positive / sum(assessed$active)
    rows[[length(rows) + 1L]] <- data.frame(
      analysis_protocol = analysis_protocol,
      permutation_rounds = rounds,
      genes = nrow(assessed),
      distinct_p_values = length(unique(p_value)),
      smallest_positive_p = min(positive),
      predicted_resolution = 1 / (nrow(assessed) * rounds),
      genes_with_zero_p = sum(p_value == 0),
      nominal_fdr = cutoff,
      calls = sum(called),
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
  }
}
resolution <- do.call(rbind, rows)
write.csv(
  resolution,
  file.path(
    "data", "derived", "crispulator_facs_mageck_permutation_resolution.csv"
  ),
  row.names = FALSE
)

cat("\nPermutation resolution, one MOI 0.20 realization, 10,000 genes:\n")
print(unique(resolution[, c(
  "permutation_rounds", "distinct_p_values", "smallest_positive_p",
  "predicted_resolution", "genes_with_zero_p"
)]), row.names = FALSE)
cat("\nCalls, realized FDP, and F1 by requested cutoff:\n")
print(resolution[, c(
  "permutation_rounds", "nominal_fdr", "calls", "realized_fdp", "f1"
)], row.names = FALSE)
