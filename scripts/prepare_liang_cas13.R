#!/usr/bin/env Rscript

# Prepare the public inputs for the Liang et al. transcriptome-scale Cas13
# fitness benchmark (Cell Genomics, 2026; BioProject PRJNA1344834).
#
# Large source workbooks remain under data/raw/ (git-ignored).  This script
# writes a compact guide library, expression table, published RRA results, and
# endpoint-run manifest that are consumed by examples/liang_cas13_benchmark.R.

options(stringsAsFactors = FALSE)

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("The `readxl` package is required.")
}

raw_dir <- file.path("data", "raw", "liang_cas13")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

s1_url <- paste0(
  "https://ars.els-cdn.com/content/image/",
  "1-s2.0-S2666979X26001151-mmc2.xlsx"
)
s2_url <- paste0(
  "https://ars.els-cdn.com/content/image/",
  "1-s2.0-S2666979X26001151-mmc3.xlsx"
)
runinfo_url <- paste0(
  "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?",
  "acc=PRJNA1344834"
)
ena_url <- paste0(
  "https://www.ebi.ac.uk/ena/portal/api/filereport?",
  "accession=PRJNA1344834&result=read_run&",
  "fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes&format=tsv"
)

download_if_missing <- function(url, destination) {
  if (!file.exists(destination) || file.info(destination)$size == 0) {
    message("Downloading ", basename(destination), " ...")
    status <- download.file(url, destination, mode = "wb", quiet = FALSE)
    if (!identical(status, 0L)) {
      stop("Download failed: ", url)
    }
  }
}

s1_path <- file.path(raw_dir, "liang_table_s1.xlsx")
s2_path <- file.path(raw_dir, "liang_table_s2.xlsx")
runinfo_path <- file.path(raw_dir, "PRJNA1344834_runinfo.csv")
ena_path <- file.path(raw_dir, "PRJNA1344834_ena_fastq.tsv")
download_if_missing(s1_url, s1_path)
download_if_missing(s2_url, s2_path)
download_if_missing(runinfo_url, runinfo_path)
download_if_missing(ena_url, ena_path)

read_sheet <- function(path, sheet) {
  as.data.frame(
    readxl::read_excel(path, sheet = sheet, skip = 2),
    check.names = FALSE
  )
}

guide <- read_sheet(s1_path, "S1B")
names(guide) <- c("sgrna", "gene", "sequence", "target_group")
guide <- guide[complete.cases(guide[, c("sgrna", "gene", "sequence")]), ]
guide$sequence <- toupper(guide$sequence)
if (nrow(guide) != 56322L || anyDuplicated(guide$sgrna)) {
  stop("Unexpected Table S1B guide library.")
}
if (any(nchar(guide$sequence) != 23L)) {
  stop("All Liang guide sequences were expected to be 23 nucleotides.")
}

write.table(
  guide,
  file.path(raw_dir, "guide_library.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
reference_path <- file.path(raw_dir, "guide_reference.fa")
reference_connection <- file(reference_path, open = "wt")
on.exit(close(reference_connection), add = TRUE)
writeLines(
  as.vector(rbind(paste0(">", guide$sgrna), guide$sequence)),
  reference_connection
)
close(reference_connection)
on.exit(NULL, add = FALSE)

# Combine total-RNA and mRNA measurements.  The paper's downstream screen
# filter is TPM > 0; using the maximum available assay avoids declaring a gene
# unexpressed merely because one library type did not detect it.
read_expression <- function(sheet) {
  x <- read_sheet(s1_path, sheet)
  names(x)[1L] <- "gene"
  for (column in setdiff(names(x), "gene")) {
    x[[column]] <- suppressWarnings(as.numeric(x[[column]]))
  }
  x
}

total_lnc <- read_expression("S1C")
mrna_lnc <- read_expression("S1E")
cell_lines <- c("HAP1", "HEK293FT", "K562", "MDA-MB-231", "THP1")

pick_expression_column <- function(x, cell_line) {
  candidates <- grep(
    paste0("^", gsub("-", "[-]", cell_line), "( RfxCas13d)?$"),
    names(x),
    value = TRUE
  )
  if (!length(candidates)) {
    return(rep(NA_real_, nrow(x)))
  }
  engineered <- grep("RfxCas13d$", candidates, value = TRUE)
  x[[if (length(engineered)) engineered[1L] else candidates[1L]]]
}

expression <- do.call(rbind, lapply(cell_lines, function(cell_line) {
  genes <- union(total_lnc$gene, mrna_lnc$gene)
  total_value <- pick_expression_column(total_lnc, cell_line)[
    match(genes, total_lnc$gene)
  ]
  mrna_value <- pick_expression_column(mrna_lnc, cell_line)[
    match(genes, mrna_lnc$gene)
  ]
  value <- pmax(total_value, mrna_value, na.rm = TRUE)
  value[!is.finite(value)] <- NA_real_
  data.frame(cell_line = cell_line, gene = genes, tpm = value)
}))
write.table(
  expression,
  file.path(raw_dir, "lncrna_expression.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

published_sheets <- setNames(
  paste0("S2", LETTERS[6:10]),
  cell_lines
)
published <- do.call(rbind, lapply(cell_lines, function(cell_line) {
  x <- read_sheet(s2_path, published_sheets[[cell_line]])
  names(x) <- c(
    "gene", "day7_p_value", "day7_log2_fold_change",
    "day14_p_value", "day14_log2_fold_change"
  )
  x$cell_line <- cell_line
  x[, c(
    "cell_line", "gene", "day7_p_value", "day7_log2_fold_change",
    "day14_p_value", "day14_log2_fold_change"
  )]
}))
write.table(
  published,
  file.path(raw_dir, "published_liang_rra.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Keep the deposited normalized counts for the primary same-input sensitivity
# analysis and for resolving the K562 day-0 sample mismatch.  The benchmark
# script rounds these fractional processed values once and gives the identical
# pseudo-count matrix to BARCS and both official MAGeCK workflows.
count_sheets <- setNames(paste0("S2", LETTERS[1:5]), cell_lines)
for (cell_line in cell_lines) {
  x <- read_sheet(s2_path, count_sheets[[cell_line]])
  names(x)[1L] <- "sgrna"
  write.table(
    x,
    file.path(
      raw_dir,
      paste0("published_processed_counts_", gsub("-", "_", cell_line), ".tsv")
    ),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

runinfo <- read.csv(runinfo_path, check.names = FALSE)
manifest <- runinfo[
  grepl("_(Day00|Day14)_R[12]$", runinfo$LibraryName),
  c("Run", "LibraryName", "spots", "size_MB", "download_path")
]
names(manifest) <- c(
  "run", "sample", "spots", "size_mb", "sra_download_path"
)
manifest$cell_line <- sub("_(Day00|Day14)_R[12]$", "", manifest$sample)
manifest$day <- as.integer(sub(
  ".*_Day([0-9]{2})_R[12]$", "\\1", manifest$sample
))
manifest$replicate <- as.integer(sub(".*_R([12])$", "\\1", manifest$sample))
ena <- read.delim(ena_path, check.names = FALSE)
ena_url_by_run <- setNames(ena$fastq_ftp, ena$run_accession)
manifest$ena_fastq_url <- paste0(
  "https://",
  ena_url_by_run[manifest$run]
)
if (any(!nzchar(manifest$ena_fastq_url)) ||
    any(is.na(manifest$ena_fastq_url))) {
  stop("ENA did not report one FASTQ URL for every endpoint run.")
}
manifest <- manifest[
  order(manifest$cell_line, manifest$day, manifest$replicate),
]
if (nrow(manifest) != 20L) {
  stop("Expected 20 day-0/day-14 Cas13 amplicon runs; found ", nrow(manifest))
}
write.table(
  manifest,
  file.path(raw_dir, "endpoint_sample_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message(
  "Prepared ", nrow(guide), " guides, ", nrow(published),
  " published gene results, and ", nrow(manifest), " endpoint runs in ",
  raw_dir, "."
)
