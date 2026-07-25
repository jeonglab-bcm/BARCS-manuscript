#!/usr/bin/env Rscript

# Generate compact Liang et al. HAP1 examples for BARCS Studio.
#
# The count example is a selected-guide view of the deposited normalized,
# ComBat-corrected values from Table S2, rounded exactly as in the benchmark.
# Its metadata retains totals calculated across the complete 56,322-guide
# table. The FASTQ files are synthetic teaching fixtures: they use the real
# selected guide sequences and emit one exact-matching read per rounded
# selected-guide pseudo-count. They are not deposited sequencing reads.

options(stringsAsFactors = FALSE)

raw_dir <- file.path("data", "raw", "liang_cas13")
output_dir <- file.path("web", "public", "examples", "liang-hap1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

guide_path <- file.path(raw_dir, "guide_library.tsv")
count_path <- file.path(raw_dir, "published_processed_counts_HAP1.tsv")
expression_path <- file.path(raw_dir, "lncrna_expression.tsv")
published_path <- file.path(raw_dir, "published_liang_rra.tsv")
required <- c(guide_path, count_path, expression_path, published_path)
if (any(!file.exists(required))) {
  stop(
    "Liang source files are missing. Run ",
    "`Rscript scripts/prepare_liang_cas13.R` first."
  )
}

guide <- read.delim(guide_path, check.names = FALSE)
processed <- read.delim(count_path, check.names = FALSE)
expression <- read.delim(expression_path, check.names = FALSE)
published <- read.delim(published_path, check.names = FALSE)

reverse_complement <- function(sequence) {
  paste(rev(chartr(
    "ACGT",
    "TGCA",
    strsplit(sequence, "", fixed = TRUE)[[1L]]
  )), collapse = "")
}

count_columns <- grep(
  "^Day(0|14) Replicate [12] +\\(Count\\)$",
  names(processed),
  value = TRUE
)
if (length(count_columns) != 4L) {
  stop("Expected four HAP1 endpoint count columns.")
}
sample_names <- c(
  "HAP1_Day00_R1", "HAP1_Day00_R2",
  "HAP1_Day14_R1", "HAP1_Day14_R2"
)

guide_index <- match(processed$sgrna, guide$sgrna)
if (anyNA(guide_index)) {
  stop("Every processed guide must occur in the Liang guide library.")
}
processed$gene <- guide$gene[guide_index]
processed$sequence <- guide$sequence[guide_index]
processed$target_group <- guide$target_group[guide_index]
canonical_sequence <- vapply(
  guide$sequence,
  function(sequence) min(sequence, reverse_complement(sequence)),
  character(1L)
)
canonical_frequency <- table(canonical_sequence)
processed_canonical <- vapply(
  processed$sequence,
  function(sequence) min(sequence, reverse_complement(sequence)),
  character(1L)
)
processed$unique_sequence <- canonical_frequency[processed_canonical] == 1L
processed$day0_total <- rowSums(processed[
  grep("^Day0.*\\(Count\\)$", names(processed))
])

target_class <- c(
  PRPF8 = "essential protein-coding",
  RPLP2 = "essential protein-coding",
  EIF3B = "essential protein-coding",
  SFPQ = "essential protein-coding",
  Hum_XLOC_053432 = "published lncRNA signal",
  Hum_XLOC_004254 = "published lncRNA signal",
  Hum_XLOC_000094 = "TPM-zero lncRNA null",
  Hum_XLOC_000105 = "TPM-zero lncRNA null"
)

select_target <- function(gene_name) {
  rows <- processed[
    processed$gene == gene_name & processed$unique_sequence,
  ]
  rows <- rows[order(-rows$day0_total, rows$sgrna), ]
  if (nrow(rows) < 4L) {
    stop(gene_name, " has fewer than four uniquely alignable guides.")
  }
  rows[seq_len(4L), ]
}

selected_targets <- do.call(
  rbind,
  lapply(names(target_class), select_target)
)
selected_controls <- processed[
  processed$target_group == "non-targeting" & processed$unique_sequence,
]
selected_controls <- selected_controls[
  order(-selected_controls$day0_total, selected_controls$sgrna),
]
if (nrow(selected_controls) < 40L) {
  stop("Fewer than 40 uniquely alignable non-targeting guides were found.")
}
selected_controls <- selected_controls[seq_len(40L), ]
selected_controls$gene <- sprintf(
  "NTC_%02d",
  rep(seq_len(10L), each = 4L)
)
selected <- rbind(selected_targets, selected_controls)
selected$control <- selected$target_group == "non-targeting"
selected$source_class <- ifelse(
  selected$control,
  "non-targeting control",
  unname(target_class[selected$gene])
)

rounded_all <- round(as.matrix(processed[count_columns]))
storage.mode(rounded_all) <- "integer"
full_library_totals <- colSums(rounded_all)
rounded_selected <- round(as.matrix(selected[count_columns]))
storage.mode(rounded_selected) <- "integer"
colnames(rounded_selected) <- sample_names

count_example <- data.frame(
  guide = selected$sgrna,
  gene = selected$gene,
  control = tolower(as.character(selected$control)),
  rounded_selected,
  check.names = FALSE
)
write.csv(
  count_example,
  file.path(output_dir, "liang-hap1-counts.csv"),
  row.names = FALSE,
  quote = FALSE
)

metadata <- data.frame(
  sample = sample_names,
  day14 = c(0L, 0L, 1L, 1L),
  replicate = c("R1", "R2", "R1", "R2"),
  total = as.integer(full_library_totals)
)
write.csv(
  metadata,
  file.path(output_dir, "liang-hap1-metadata.csv"),
  row.names = FALSE,
  quote = FALSE
)

library_example <- data.frame(
  guide = selected$sgrna,
  gene = selected$gene,
  sequence = selected$sequence,
  control = tolower(as.character(selected$control)),
  source_class = selected$source_class,
  original_target_group = selected$target_group
)
write.csv(
  library_example,
  file.path(output_dir, "liang-hap1-guide-library.csv"),
  row.names = FALSE,
  quote = FALSE
)

fastq_metadata <- metadata[, c("sample", "day14", "replicate")]
write.csv(
  fastq_metadata,
  file.path(output_dir, "liang-hap1-fastq-metadata.csv"),
  row.names = FALSE,
  quote = FALSE
)
write.csv(
  count_example,
  file.path(output_dir, "liang-hap1-fastq-expected-counts.csv"),
  row.names = FALSE,
  quote = FALSE
)

for (column in seq_along(sample_names)) {
  path <- file.path(output_dir, paste0(sample_names[column], ".fastq.gz"))
  connection <- gzfile(path, open = "wt", compression = 9)
  read_number <- 0L
  for (row in seq_len(nrow(selected))) {
    count <- rounded_selected[row, column]
    if (!count) next
    spacer <- selected$sequence[row]
    for (copy in seq_len(count)) {
      read_number <- read_number + 1L
      reverse <- read_number %% 10L == 0L
      read_sequence <- if (reverse) {
        reverse_complement(spacer)
      } else {
        spacer
      }
      writeLines(
        c(
          sprintf("@%s_%07d", sample_names[column], read_number),
          read_sequence,
          "+",
          paste(rep("I", nchar(read_sequence)), collapse = "")
        ),
        connection
      )
    }
  }
  close(connection)
}

hap1_expression <- expression[expression$cell_line == "HAP1", ]
hap1_published <- published[published$cell_line == "HAP1", ]
provenance <- unique(selected[, c("gene", "source_class")])
provenance$tpm <- hap1_expression$tpm[
  match(provenance$gene, hap1_expression$gene)
]
provenance$published_day14_p_value <- hap1_published$day14_p_value[
  match(provenance$gene, hap1_published$gene)
]
provenance$published_day14_log2_fold_change <-
  hap1_published$day14_log2_fold_change[
    match(provenance$gene, hap1_published$gene)
  ]
write.csv(
  provenance,
  file.path(output_dir, "liang-hap1-provenance.csv"),
  row.names = FALSE,
  quote = FALSE,
  na = ""
)

readme <- c(
  "# Liang HAP1 BARCS examples",
  "",
  "Source study: Liang et al. transcriptome-scale Cas13 fitness screens ",
  "(Cell Genomics, 2026; BioProject PRJNA1344834).",
  "",
  "## Count-matrix example",
  "",
  "`liang-hap1-counts.csv` contains 72 real Liang guide rows: four guides ",
  "each for four essential/protein-coding genes, two published lncRNA ",
  "signals, two TPM-zero lncRNA nulls, and 40 non-targeting controls. Values ",
  "are Liang's deposited normalized, ComBat-corrected processed counts rounded ",
  "to integer pseudo-counts. `liang-hap1-metadata.csv` preserves full-library ",
  "totals calculated before the 72-guide subset was selected.",
  "",
  "This is a compact same-input sensitivity example. It is not a raw-count ",
  "likelihood analysis and should not be used to replace the complete benchmark.",
  "",
  "## FASTQ example",
  "",
  "The four `.fastq.gz` files are synthetic teaching fixtures, not Liang's ",
  "deposited sequencing reads. They use the real selected 23-nt Liang guide ",
  "sequences and contain one exact-matching read per rounded selected-guide ",
  "pseudo-count. Ten percent of reads contain the reverse-complement spacer to ",
  "exercise orientation detection. Alignment should reproduce ",
  "`liang-hap1-fastq-expected-counts.csv` exactly.",
  "",
  "Use `liang-hap1-guide-library.csv` as the guide library and ",
  "`liang-hap1-fastq-metadata.csv` as sample metadata. Fit `day14` while ",
  "adjusting for `replicate`.",
  "",
  "The FASTQ-derived analysis uses totals from the selected synthetic library; ",
  "therefore its estimates need not equal the count example, which retains the ",
  "complete Liang library denominators."
)
writeLines(readme, file.path(output_dir, "README.md"))

bundle_files <- c(
  "README.md",
  "liang-hap1-guide-library.csv",
  "liang-hap1-fastq-metadata.csv",
  "liang-hap1-fastq-expected-counts.csv",
  paste0(sample_names, ".fastq.gz")
)
old_directory <- setwd(output_dir)
on.exit(setwd(old_directory), add = TRUE)
unlink("liang-hap1-fastq-example.zip")
utils::zip(
  "liang-hap1-fastq-example.zip",
  bundle_files,
  flags = "-q -j"
)
setwd(old_directory)
on.exit(NULL, add = FALSE)

message(
  "Wrote ", nrow(selected), " guides and ",
  sum(rounded_selected), " synthetic reads to ", output_dir, "."
)
