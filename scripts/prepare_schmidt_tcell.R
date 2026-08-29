#!/usr/bin/env Rscript

# Prepare the public model inputs for the Schmidt et al. primary human T cell
# CRISPRa/CRISPRi screen re-analysis (Science 2022, doi:10.1126/science.abj4008).
#
# Downloads the GEO raw sgRNA read-count workbook (GSE174255) and writes one
# counts/metadata pair per library set into results/schmidt_tcell/input/:
#   {CRISPRa_SetA,CRISPRa_SetB,CRISPRi_SetA,CRISPRi_SetB}_{counts,metadata}.csv
#
# Unlike the Liu re-analysis, these are true raw integer read counts: GEO
# publishes the per-library counts the authors aligned with MAGeCK, before any
# normalisation. The four sheets are four independently sequenced library
# pools, so each is prepared, and later fitted, on its own denominators. The
# 22 MB workbook and the ~15 MB inputs both stay under git-ignored paths; only
# the compact gene-level results are committed.

options(stringsAsFactors = FALSE)

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("The `readxl` package is required.")
}

raw_dir <- file.path("data", "raw", "schmidt_tcell")
input_dir <- file.path("results", "schmidt_tcell", "input")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

counts_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174255/suppl/",
  "GSE174255_sgRNA-Read-Counts.xlsx"
)
workbook_path <- file.path(raw_dir, "GSE174255_sgRNA-Read-Counts.xlsx")

if (!file.exists(workbook_path) || file.info(workbook_path)$size == 0) {
  message("Downloading GSE174255 sgRNA read counts ...")
  if (!identical(download.file(counts_url, workbook_path, mode = "wb"), 0L)) {
    stop("Download failed: ", counts_url)
  }
}

# The workbook prefixes every sheet with "CRISPRa", but Calabrese is the CRISPRa
# library and Dolcetto the CRISPRi library (Methods, and the GEO sample titles).
# Modality is taken from the library, not from the sheet name.
library_sets <- list(
  CRISPRa_SetA = list(sheet = "CRISPRa.CalabreseSetA.count", modality = "CRISPRa"),
  CRISPRa_SetB = list(sheet = "CRISPRa.CalabreseSetB.count", modality = "CRISPRa"),
  CRISPRi_SetA = list(sheet = "CRISPRa.DolcettoSetA.count",  modality = "CRISPRi"),
  CRISPRi_SetB = list(sheet = "CRISPRa.DolcettoSetB.count",  modality = "CRISPRi")
)

for (set_name in names(library_sets)) {
  spec <- library_sets[[set_name]]
  sheet <- as.data.frame(
    suppressMessages(readxl::read_excel(workbook_path, sheet = spec$sheet)),
    check.names = FALSE
  )

  # The plasmid library is the pre-transduction pool, not a sorted phenotype
  # bin; it plays no part in the high-versus-low contrast and is dropped.
  sample_columns <- setdiff(colnames(sheet), c("sgRNA", "Gene"))
  sample_columns <- grep("_Plasmid$", sample_columns, value = TRUE, invert = TRUE)

  count_matrix <- as.matrix(sheet[, sample_columns, drop = FALSE])
  storage.mode(count_matrix) <- "double"
  if (any(abs(count_matrix - round(count_matrix)) > 0) || any(count_matrix < 0)) {
    stop("GEO counts must be non-negative integers; sheet ", spec$sheet)
  }

  counts <- data.frame(
    guide = sheet$sgRNA,
    gene = sheet$Gene,
    # The library's 496 nontargeting sgRNAs, the controls the published MAGeCK
    # analysis passed to --control-sgrna.
    control = ifelse(sheet$Gene == "NO-TARGET", "true", "false"),
    count_matrix,
    check.names = FALSE
  )

  metadata <- data.frame(
    sample = sample_columns,
    library_set = set_name,
    modality = spec$modality,
    donor = sub(".*_(Donor[12])_.*", "\\1", sample_columns),
    assay = sub(".*_(IFNG|IL2)_.*", "\\1", sample_columns),
    bin = sub(".*_(low|high|unsorted)$", "\\1", sample_columns)
  )
  if (anyNA(metadata) || !all(metadata$bin %in% c("low", "high", "unsorted"))) {
    stop("Could not decode every sample name in sheet ", spec$sheet)
  }
  # Beta-binomial denominators: reads sequenced from this library, which is the
  # column sum within its own pool. No cross-set rescaling is applied.
  metadata$total <- colSums(count_matrix)[metadata$sample]

  write.csv(
    counts, file.path(input_dir, paste0(set_name, "_counts.csv")),
    row.names = FALSE, quote = FALSE
  )
  write.csv(
    metadata, file.path(input_dir, paste0(set_name, "_metadata.csv")),
    row.names = FALSE, quote = FALSE
  )
  message(sprintf(
    "%s: %d guides (%d nontargeting) x %d libraries, %.1fM reads",
    set_name, nrow(counts), sum(counts$control == "true"),
    length(sample_columns), sum(count_matrix) / 1e6
  ))
}
