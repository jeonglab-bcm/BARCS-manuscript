#!/usr/bin/env Rscript

# Generate full and compact Liang et al. HAP1 examples for BARCS Studio.
#
# The primary count example contains every HAP1 guide present in the deposited
# normalized, ComBat-corrected Table S2 values, rounded exactly as in the
# benchmark. The small FASTQ files remain explicit software-validation
# fixtures. A separate manifest points to the four deposited HAP1 endpoint
# FASTQs; no synthetic reads are presented as biological data.

options(stringsAsFactors = FALSE)

raw_dir <- file.path("data", "raw", "liang_cas13")
output_dir <- file.path("web", "public", "examples", "liang-hap1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

guide_path <- file.path(raw_dir, "guide_library.tsv")
count_path <- file.path(raw_dir, "published_processed_counts_HAP1.tsv")
expression_path <- file.path(raw_dir, "lncrna_expression.tsv")
published_path <- file.path(raw_dir, "published_liang_rra.tsv")
endpoint_manifest_path <- file.path(raw_dir, "endpoint_sample_manifest.tsv")
ena_manifest_path <- file.path(raw_dir, "PRJNA1344834_ena_fastq.tsv")
required <- c(
  guide_path, count_path, expression_path, published_path,
  endpoint_manifest_path, ena_manifest_path
)
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
endpoint_manifest <- read.delim(endpoint_manifest_path, check.names = FALSE)
ena_manifest <- read.delim(ena_manifest_path, check.names = FALSE)

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

# The complete browser example uses every guide for which Liang deposited HAP1
# processed endpoint values. Table S2 omits 148 of the 56,322 library entries,
# so this file has 56,174 rows rather than silently filling absent values.
full_gene <- processed$gene
full_control <- processed$target_group == "non-targeting"
control_index <- which(full_control)
full_gene[control_index] <- sprintf(
  "NTC_%03d",
  ceiling(seq_along(control_index) / 4)
)
full_count_example <- data.frame(
  guide = processed$sgrna,
  gene = full_gene,
  control = tolower(as.character(full_control)),
  rounded_all,
  check.names = FALSE
)
colnames(full_count_example)[-(1:3)] <- sample_names
write.csv(
  full_count_example,
  file.path(output_dir, "liang-hap1-full-counts.csv"),
  row.names = FALSE,
  quote = FALSE
)

raw_library_control <- guide$target_group == "non-targeting"
raw_library_gene <- guide$gene
raw_control_index <- which(raw_library_control)
raw_library_gene[raw_control_index] <- sprintf(
  "NTC_%03d",
  ceiling(seq_along(raw_control_index) / 4)
)
full_library_example <- data.frame(
  guide = guide$sgrna,
  gene = raw_library_gene,
  sequence = guide$sequence,
  control = tolower(as.character(raw_library_control)),
  original_gene = guide$gene,
  target_group = guide$target_group
)
write.csv(
  full_library_example,
  file.path(output_dir, "liang-hap1-full-guide-library.csv"),
  row.names = FALSE,
  quote = FALSE
)

hapi_endpoint <- endpoint_manifest[
  endpoint_manifest$cell_line == "HAP1" &
    endpoint_manifest$day %in% c(0L, 14L),
]
hapi_endpoint <- hapi_endpoint[
  match(sample_names, hapi_endpoint$sample),
]
ena_index <- match(hapi_endpoint$run, ena_manifest$run_accession)
if (nrow(hapi_endpoint) != 4L || anyNA(ena_index)) {
  stop("Could not resolve the four deposited HAP1 endpoint FASTQs.")
}
real_fastq_manifest <- data.frame(
  sample = hapi_endpoint$sample,
  run_accession = hapi_endpoint$run,
  day14 = as.integer(hapi_endpoint$day == 14L),
  replicate = paste0("R", hapi_endpoint$replicate),
  downloaded_filename = paste0(hapi_endpoint$run, ".fastq.gz"),
  browser_filename = paste0(hapi_endpoint$sample, ".fastq.gz"),
  fastq_url = hapi_endpoint$ena_fastq_url,
  md5 = ena_manifest$fastq_md5[ena_index],
  compressed_bytes = as.numeric(ena_manifest$fastq_bytes[ena_index])
)
write.csv(
  real_fastq_manifest,
  file.path(output_dir, "liang-hap1-real-fastq-manifest.csv"),
  row.names = FALSE,
  quote = FALSE
)
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
  "# Liang HAP1 full real-data case",
  "",
  "Source study: Liang et al. transcriptome-scale Cas13 fitness screens ",
  "(Cell Genomics, 2026; BioProject PRJNA1344834).",
  "",
  "## Full processed-count analysis",
  "",
  "`liang-hap1-full-counts.csv` contains all 56,174 HAP1 guides for which ",
  "Liang deposited processed endpoint values. The four columns are the two ",
  "day-0 and two day-14 libraries. Values are the deposited normalized, ",
  "ComBat-corrected values rounded to integer pseudo-counts, exactly as in ",
  "the manuscript sensitivity analysis. `liang-hap1-metadata.csv` retains ",
  "the corresponding complete-library totals. Non-targeting guides are ",
  "assigned to four-guide pseudo-genes only for gene-level null calibration; ",
  "their guide identities and counts are unchanged.",
  "",
  "This is a complete real screen, not a selected-gene toy case. Because the ",
  "deposited values have already been normalized and corrected, it remains a ",
  "processed-count sensitivity analysis rather than a raw-count likelihood ",
  "analysis.",
  "",
  "## Deposited FASTQ analysis",
  "",
  "`liang-hap1-real-fastq-manifest.csv` lists the four real HAP1 endpoint ",
  "FASTQs from ENA, including run accessions, HTTPS URLs, MD5 checksums, and ",
  "compressed byte sizes. Together they are about 1.7 GB. ",
  "`liang-hap1-full-guide-library.csv` contains all 56,322 library guides. ",
  "For browser analysis, rename each downloaded file to the ",
  "`browser_filename` in the manifest so it matches the sample metadata. ",
  "From the repository root, run:",
  "",
  "```sh",
  "bash scripts/run_liang_hap1_real_case.sh",
  "```",
  "",
  "The runner streams the deposited reads through the article-matched ",
  "Cutadapt/Bowtie workflow and writes a BARCS-ready raw count matrix. It does ",
  "not retain the FASTQs.",
  "",
  "## Small validation fixture",
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
  "The four small `.fastq.gz` files are synthetic software-test fixtures, not ",
  "a biological example and not Liang's ",
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
writeLines(
  trimws(readme, which = "right"),
  file.path(output_dir, "README.md")
)

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
  "Wrote the ", nrow(full_count_example), "-guide real HAP1 case plus ",
  nrow(selected), " validation guides and ",
  sum(rounded_selected), " synthetic test reads to ", output_dir, "."
)
