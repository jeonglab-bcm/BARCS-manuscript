#!/usr/bin/env Rscript

# MAUDE calibration on GSE242880 held-out non-targeting guides.
#
# MAUDE consumes non-targeting guides to build its empirical null, so its
# control calibration cannot be measured on the same guides it was given: those
# Z scores are centred by construction. We therefore split the 593 non-targeting
# guides in half. One half is declared to MAUDE as `negativeControl`; the other
# half is relabelled so each held-out guide becomes its own element and is
# scored like any other target. Because a held-out non-targeting guide has no
# true effect, the rate at which it is called is a direct measurement of
# MAUDE's realised error rate against its nominal one.
#
# MAUDE also requires an unsorted reference sample (`unsortedBin`). The
# low-coverage, high-MOI arm used for the BARCS benchmark has only Q1--Q4 and no
# unsorted column, so this analysis uses the high-coverage, high-MOI arm, which
# deposits a GFP unsorted sample for each donor. That is a different arm from
# the BARCS comparison and the two numbers are not interchangeable.
#
#     Rscript examples/waterbear_maude_calibration.R

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(MAUDE))

count_path <- file.path(
  "data", "raw", "waterbear",
  "GSE242880_high_low_MOI_screens_D1_D2_D3_2022-04-23.count.txt.gz"
)
if (!file.exists(count_path)) {
  stop("GSE242880 counts are missing; see data/raw/waterbear/README.md.", call. = FALSE)
}

raw <- read.delim(gzfile(count_path), check.names = FALSE)
donors <- c("D1", "D2", "D3")
bin_masses <- c(0.2, 0.3, 0.3, 0.2)
control_label <- "Non-Targeting Control"
seed <- 20250724

counts <- do.call(rbind, lapply(donors, function(donor) {
  sorted <- sprintf("High_Coverage_High_MOI_%s_Q%d", donor, 1:4)
  unsorted <- sprintf("High_Coverage_High_MOI_%s_GFP", donor)
  if (!all(c(sorted, unsorted) %in% names(raw))) {
    stop("Expected four sorted bins and one GFP column for ", donor, ".")
  }
  data.frame(
    sgRNA = raw$sgRNA,
    element = raw$Gene,
    NT = raw$Gene == control_label,
    A = raw[[sorted[1]]], B = raw[[sorted[2]]],
    C = raw[[sorted[3]]], D = raw[[sorted[4]]],
    NS = raw[[unsorted]],
    expt = donor
  )
}))

bin_model <- MAUDE::makeBinModel(
  data.frame(Bin = c("A", "B", "C", "D"), fraction = bin_masses)
)
bin_model <- bin_model[bin_model$Bin %in% c("A", "B", "C", "D"), ]
bin_stats <- do.call(rbind, lapply(donors, function(donor) {
  stats <- bin_model
  stats$expt <- donor
  stats
}))

set.seed(seed)
control_guides <- sort(unique(counts$sgRNA[counts$NT]))
calibrating <- sample(control_guides, floor(length(control_guides) / 2))
held_out <- setdiff(control_guides, calibrating)

split_counts <- counts
split_counts$NT <- split_counts$sgRNA %in% calibrating
split_counts$element <- ifelse(
  split_counts$sgRNA %in% held_out,
  paste0("HELDOUT_", split_counts$sgRNA),
  split_counts$element
)

guide_stats <- MAUDE::findGuideHitsAllScreens(
  experiments = unique(split_counts["expt"]),
  countDataFrame = split_counts,
  binStats = bin_stats,
  sortBins = c("A", "B", "C", "D"),
  unsortedBin = "NS",
  negativeControl = "NT"
)
element_stats <- MAUDE::getElementwiseStats(
  experiments = unique(guide_stats["expt"]),
  normNBSummaries = guide_stats,
  elementIDs = "element",
  negativeControl = "NT"
)
element_stats$fdr <- p.adjust(element_stats$p.value, method = "BH")

is_held <- grepl("^HELDOUT_", element_stats$element)
held_stats <- element_stats[is_held, ]
target_stats <- element_stats[!is_held, ]

thresholds <- c(0.01, 0.05, 0.10, 0.20)
calibration <- data.frame(
  arm = "High_Coverage_High_MOI",
  method = "MAUDE",
  control_guides_total = length(control_guides),
  control_guides_calibrating = length(calibrating),
  control_guides_held_out = length(held_out),
  held_out_tests = nrow(held_stats),
  nominal_threshold = thresholds,
  observed_rate = sapply(thresholds, function(t) mean(held_stats$p.value < t)),
  ratio_to_nominal = sapply(
    thresholds, function(t) mean(held_stats$p.value < t) / t
  ),
  seed = seed
)
summary_row <- data.frame(
  arm = "High_Coverage_High_MOI",
  method = "MAUDE",
  held_out_called_fdr_0_10 = mean(held_stats$fdr < 0.10),
  target_elements = nrow(target_stats),
  target_called_fdr_0_10 = sum(target_stats$fdr < 0.10),
  seed = seed
)

dir.create(file.path("data", "derived"), showWarnings = FALSE, recursive = TRUE)
write.csv(
  calibration,
  file.path("data", "derived", "waterbear_maude_heldout_calibration.csv"),
  row.names = FALSE
)
write.csv(
  summary_row,
  file.path("data", "derived", "waterbear_maude_heldout_summary.csv"),
  row.names = FALSE
)

print(calibration)
print(summary_row)
