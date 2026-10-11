#!/usr/bin/env Rscript

# Prepare the genome-wide model inputs for the Liu et al. in vivo T cell
# IFNγ screen re-analysis (Nature 2026, doi:10.1038/s41586-026-10906-9).
#
# Supplementary Table 2 is the genome-wide IFNγ screen behind the paper's
# main Fig. 4 (19,113 genes; MAGeCK RRA on IFNγ-high vs IFNγ-low sorted
# tumour-infiltrating T cells, n = 2 donors). Its sgRNA-summary sheet stores
# one row per guide-and-donor: the `control_count` field is the IFNγ-low gate
# and `treatment_count` the IFNγ-high gate, and the `_r0` / `_r1` suffix on the
# sgRNA identifier is the donor replicate. This script reshapes that long
# table into the guide-by-library count matrix and metadata BARCS needs, with
# a `~ gate + donor` design (four libraries per guide: two gates x two donors).
#
# Writes:
#   results/liu_genomewide/input/counts.csv
#   results/liu_genomewide/input/metadata.csv
#
# Unlike the focused sub-library (Supplementary Table 5, scripts/prepare_liu_tcell.R),
# the genome-wide screen is collapsed to two donor replicates, so only a
# donor-adjusted fit is possible here; there is no per-mouse structure to pair
# on. The large workbook stays under data/raw/ (git-ignored); only the parsed
# inputs are committed.

options(stringsAsFactors = FALSE)

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("The `readxl` package is required.")
}

raw_dir <- file.path("data", "raw", "liu_tcell")
input_dir <- file.path("results", "liu_genomewide", "input")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

table2_path <- file.path(raw_dir, "Supplementary_Table_2.xlsx")
if (!file.exists(table2_path)) {
  stop(
    "Missing ", table2_path,
    ". Run `Rscript scripts/prepare_liu_tcell.R` first to fetch the archive, ",
    "or download the Nature supplementary archive into ", raw_dir, ".",
    call. = FALSE
  )
}

# MAGeCK median normalisation yields one-decimal counts; carry them to integers
# with round-half-up (base R's round() is round-half-to-even on .5 ties).
round_half_up <- function(x) floor(x + 0.5)

sgrna <- as.data.frame(
  suppressMessages(readxl::read_excel(
    table2_path, sheet = "IFNghi vs. IFNglo_sgRNA_summary"
  )),
  check.names = FALSE
)

# The donor replicate is the trailing _r0 / _r1 on the sgRNA identifier; the
# guide identity is what remains once it is stripped.
sgrna$donor <- sub(".*_(r[0-9]+)$", "\\1", sgrna$sgrna)
sgrna$guide <- sub("_r[0-9]+$", "", sgrna$sgrna)
donors <- sort(unique(sgrna$donor))
if (!all(donors == c("r0", "r1"))) {
  stop("Expected exactly two donor replicates r0 and r1.")
}

# Reshape long (guide x donor rows) to wide (one row per guide). Each guide
# must appear once per donor; the low/high gates become four named libraries.
guides <- sort(unique(sgrna$guide))
gate_field <- c(lo = "control_count", hi = "treatment_count")
library_names <- as.vector(t(outer(
  donors, names(gate_field), function(d, g) paste0(d, "_", g)
)))

count_matrix <- matrix(
  NA_real_, nrow = length(guides), ncol = length(library_names),
  dimnames = list(guides, library_names)
)
for (donor in donors) {
  block <- sgrna[sgrna$donor == donor, ]
  index <- match(guides, block$guide)
  if (anyNA(index)) {
    stop("A guide is missing donor ", donor, " in the sgRNA summary.")
  }
  for (gate in names(gate_field)) {
    count_matrix[, paste0(donor, "_", gate)] <-
      round_half_up(block[[gate_field[[gate]]]][index])
  }
}
if (anyNA(count_matrix) || any(count_matrix < 0)) {
  stop("Reshaped counts must be finite and non-negative.")
}

# One gene label per guide; NTCTRL marks the 1,000 non-targeting controls.
gene <- sgrna$Gene[match(guides, sgrna$guide)]
counts <- data.frame(
  guide = guides,
  gene = gene,
  control = ifelse(gene == "NTCTRL", "true", "false"),
  count_matrix,
  check.names = FALSE,
  row.names = NULL
)

metadata <- data.frame(
  sample = library_names,
  gate = ifelse(grepl("_hi$", library_names), 1L, 0L),
  donor = sub("_(lo|hi)$", "", library_names)
)
# Library totals are the column sums of the normalised counts, not true
# sequencing depths; they are the beta-binomial denominators.
metadata$total <- colSums(count_matrix)[metadata$sample]
metadata <- metadata[, c("sample", "gate", "donor", "total")]

write.csv(
  counts, file.path(input_dir, "counts.csv"), row.names = FALSE, quote = FALSE
)
write.csv(
  metadata, file.path(input_dir, "metadata.csv"),
  row.names = FALSE, quote = FALSE
)
message(sprintf(
  "genome-wide: %d guides x %d libraries, %d genes (%d NTCTRL guides)",
  nrow(counts), length(library_names), length(unique(gene)),
  sum(gene == "NTCTRL")
))
