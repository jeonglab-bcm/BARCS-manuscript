#!/usr/bin/env Rscript

# Prepare the public model inputs for the Liu et al. in vivo T cell CRISPR
# screen re-analysis (Nature 2026, doi:10.1038/s41586-026-10906-9).
#
# Downloads the Nature supplementary archive, reads Supplementary Table 5, and
# writes the four input files consumed by examples/liu_tcell_barcs.R:
#   results/liu_tcell/input/{armA,armB}_{counts,metadata}.csv
#
# Supplementary Table 5 is the only screen in the paper reporting per-mouse
# counts. Its per-guide `control_count` / `treatment_count` fields hold the
# MAGeCK median-normalised counts for the IFNγ-low and IFNγ-high sort gates,
# one `/`-separated value per library. Arm A (CD3-scFv vs A375low) comes from
# the pooled three-donor sheet and Arm B (NY-ESO-1 TCR vs WT A375) from the
# pooled two-donor sheet; the per-donor and per-mouse sheets are separately
# renormalised and are not used. The large source workbook stays under
# data/raw/ (git-ignored); only the compact parsed inputs are committed.

options(stringsAsFactors = FALSE)

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("The `readxl` package is required.")
}

raw_dir <- file.path("data", "raw", "liu_tcell")
input_dir <- file.path("results", "liu_tcell", "input")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

supplement_url <- paste0(
  "https://media.springernature.com/original/springer-static/esm/",
  "art%3A10.1038%2Fs41586-026-10906-9/MediaObjects/",
  "41586_2026_10906_MOESM3_ESM.zip"
)
archive_path <- file.path(raw_dir, "41586_2026_10906_MOESM3_ESM.zip")
table5_path <- file.path(raw_dir, "Supplementary_Table_5.xlsx")

if (!file.exists(archive_path) || file.info(archive_path)$size == 0) {
  message("Downloading Nature supplementary archive ...")
  if (!identical(download.file(supplement_url, archive_path, mode = "wb"), 0L)) {
    stop("Download failed: ", supplement_url)
  }
}

if (!file.exists(table5_path)) {
  members <- utils::unzip(archive_path, list = TRUE)$Name
  table5_member <- grep("Supplementary_Table_5\\.xlsx$", members, value = TRUE)
  if (length(table5_member) != 1L) {
    stop("Could not locate Supplementary_Table_5.xlsx in the archive.")
  }
  extracted <- utils::unzip(archive_path, files = table5_member, exdir = raw_dir)
  file.rename(extracted, table5_path)
}

# MAGeCK reports the median-normalised counts to one decimal; the published
# values were carried to integers with round-half-up, which base R's round()
# (round-half-to-even) does not reproduce on the exact .5 ties.
round_half_up <- function(x) floor(x + 0.5)

# One arm per pooled MAGeCK run. `mice` is the per-library order inside each
# `/`-separated count field, decoded as donor (D#/W#) and mouse (…M#).
arms <- list(
  armA = list(
    sheet = "CD3scFv_A375low_3Donor_sgRNASum",
    mice = c("D1M1", "D1M2", "D1M3", "D2M1", "D2M2", "D2M3", "D3M1", "D3M2")
  ),
  armB = list(
    sheet = "WT_A375_2donors_sgRNASum",
    mice = c("W1M1", "W1M2", "W2M1", "W2M2")
  )
)

split_counts <- function(field, n_libraries) {
  values <- lapply(strsplit(field, "/", fixed = TRUE), as.numeric)
  if (any(lengths(values) != n_libraries)) {
    stop("A guide reports the wrong number of per-library counts.")
  }
  round_half_up(do.call(rbind, values))
}

for (arm_name in names(arms)) {
  arm <- arms[[arm_name]]
  sheet <- as.data.frame(
    suppressMessages(readxl::read_excel(table5_path, sheet = arm$sheet)),
    check.names = FALSE
  )

  n_libraries <- length(arm$mice)
  # control_count is the IFNγ-low gate, treatment_count the IFNγ-high gate.
  low <- split_counts(sheet$control_count, n_libraries)
  high <- split_counts(sheet$treatment_count, n_libraries)
  colnames(low) <- paste0(arm$mice, "_lo")
  colnames(high) <- paste0(arm$mice, "_hi")

  # Interleave the two gates per mouse so paired libraries sit side by side.
  sample_order <- as.vector(rbind(colnames(low), colnames(high)))
  count_matrix <- cbind(low, high)[, sample_order, drop = FALSE]

  counts <- data.frame(
    guide = sheet$sgrna,
    gene = sheet$Gene,
    control = ifelse(sheet$Gene == "NTCTRL", "true", "false"),
    count_matrix,
    check.names = FALSE
  )

  metadata <- data.frame(
    sample = sample_order,
    gate = ifelse(grepl("_hi$", sample_order), 1L, 0L),
    mouse = sub("_(lo|hi)$", "", sample_order)
  )
  # Donor is the leading platform/mouse index (D1M2 -> D1, W2M1 -> D2).
  metadata$donor <- paste0("D", substr(metadata$mouse, 2L, 2L))
  # Library totals are the column sums of the normalised counts, not true
  # sequencing depths; they are the beta-binomial denominators.
  metadata$total <- colSums(count_matrix)[metadata$sample]
  metadata <- metadata[, c("sample", "gate", "donor", "mouse", "total")]

  write.csv(
    counts, file.path(input_dir, paste0(arm_name, "_counts.csv")),
    row.names = FALSE, quote = FALSE
  )
  write.csv(
    metadata, file.path(input_dir, paste0(arm_name, "_metadata.csv")),
    row.names = FALSE, quote = FALSE
  )
  message(sprintf(
    "%s: %d guides x %d libraries -> %s_counts.csv, %s_metadata.csv",
    arm_name, nrow(counts), n_libraries, arm_name, arm_name
  ))
}
