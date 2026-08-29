# Shared pieces of the Schmidt et al. T cell screen re-analysis
# (Science 2022, doi:10.1126/science.abj4008; raw counts GEO GSE174255).
#
# Sourced by examples/schmidt_tcell_barcs.R, examples/schmidt_tcell_mageck.R,
# and examples/schmidt_tcell_null_calibration.R so that the library sets, the
# nontargeting pseudo-gene construction, and the held-out control split are
# defined once and are identical for both methods.

schmidt_screens <- c("CRISPRa_IL2", "CRISPRa_IFNG", "CRISPRi_IL2", "CRISPRi_IFNG")
schmidt_library_sets <- c("CRISPRa_SetA", "CRISPRa_SetB", "CRISPRi_SetA", "CRISPRi_SetB")

# Each library set carries its own 496 nontargeting sgRNAs. The two sets share
# no sgRNA sequence at all -- not for targeting guides and not for controls --
# so a screen has 992 distinct control guides, and the published analysis labels
# every one of them as the single gene NO-TARGET.
#
# One 992-guide gene is not a useful null: it cannot say how often a *gene-sized*
# group of null guides gets called, and it leaves `bb_gene_consistency()` with a
# single control gene, too few to estimate a control-based null at all. The
# controls are therefore grouped into pseudo-genes with the same shape as a real
# gene. Grouping runs inside each library set, three sgRNAs at a time in
# lexicographic order, and pseudo-gene k takes its three Set A guides and its
# three Set B guides -- six guides drawn from two library sets, exactly the
# structure of a real gene here. The grouping is deterministic and is shared by
# every method and every run. The final group absorbs the remainder rather than
# being left short of `min_guides`.
schmidt_pseudogene_index <- function(sgrna, group_size = 3L) {
  order_index <- match(sgrna, sort(unique(sgrna))) - 1L
  n_groups <- length(unique(sgrna)) %/% group_size
  pmin(order_index %/% group_size, n_groups - 1L)
}

# The same grouping applied to guides that have already been stacked across the
# two library sets, where the |SetA / |SetB qualifier records which set a guide
# came from.
schmidt_control_pseudogene <- function(guide) {
  set <- sub("^.*\\|", "", guide)
  index <- integer(length(guide))
  for (one_set in unique(set)) {
    in_set <- set == one_set
    index[in_set] <- schmidt_pseudogene_index(guide[in_set])
  }
  index
}

# Held-out split for the null-calibration check, assigned by pseudo-gene so that
# no pseudo-gene ever mixes calibration and evaluation guides. Alternating
# pseudo-genes go to each half, and the two folds swap the roles, so across both
# folds every control pseudo-gene is evaluated exactly once while being held out
# of the null the method was tuned on.
#
# The calibration half is the only control set a method may tune its null on --
# BARCS's control calibration, MAGeCK's --control-sgrna. The evaluation half is
# relabelled and scored alongside the real genes. Fitting a scale and evaluating
# it on the same controls is not treated as validation.
schmidt_control_fold <- function(index, fold) {
  ifelse(index %% 2L == as.integer(fold), "evaluation", "calibration")
}

schmidt_read_set <- function(input_dir, set_name) {
  counts_path <- file.path(input_dir, paste0(set_name, "_counts.csv"))
  metadata_path <- file.path(input_dir, paste0(set_name, "_metadata.csv"))
  if (!file.exists(counts_path) || !file.exists(metadata_path)) {
    stop(
      "Missing Schmidt T cell inputs: ", counts_path, " / ", metadata_path,
      "\nRun `Rscript scripts/prepare_schmidt_tcell.R` first.",
      call. = FALSE
    )
  }

  counts_table <- read.csv(counts_path, check.names = FALSE)
  metadata <- read.csv(metadata_path, check.names = FALSE)
  metadata$donor <- factor(metadata$donor)
  metadata$bin <- factor(metadata$bin, levels = c("low", "unsorted", "high"))

  count_matrix <- as.matrix(counts_table[, metadata$sample, drop = FALSE])
  storage.mode(count_matrix) <- "double"

  list(
    counts = count_matrix,
    metadata = metadata,
    guide = counts_table$guide,
    gene = counts_table$gene,
    control = tolower(as.character(counts_table$control)) == "true",
    totals = metadata$total,
    set_name = set_name
  )
}

# The published MAGeCK input for one screen, built exactly as the Methods
# describe: "raw read counts across both library sets were normalized to the
# total read count in each sample, and each of the matching samples across two
# sets were merged to generate a single normalized read count table".
#
# "Normalized to the total read count" is applied as MAGeCK's own total-count
# normalisation -- count / library total x mean library total -- so the merged
# table stays on a read-count scale, which is what `mageck test --norm-method
# none` then expects. Only the sorted bins enter the comparison; the unsorted
# libraries were not used by the published `mageck test` call.
schmidt_mageck_table <- function(input_dir, modality, assay) {
  blocks <- lapply(c("SetA", "SetB"), function(set_tag) {
    set <- schmidt_read_set(input_dir, paste(modality, set_tag, sep = "_"))
    keep <- set$metadata$assay == assay & set$metadata$bin %in% c("low", "high")

    count_matrix <- set$counts[, keep, drop = FALSE]
    totals <- colSums(count_matrix)
    normalised <- round(sweep(count_matrix, 2L, totals, "/") * mean(totals))
    # Set-independent column names, so the two blocks stack into single columns.
    colnames(normalised) <- paste(
      set$metadata$donor[keep], set$metadata$bin[keep], sep = "_"
    )

    data.frame(
      # The two sets share no sgRNA id -- the nontargeting sequences were
      # checked, not assumed -- so the merge needs no disambiguation. The set
      # is still recorded in the id, because it is what the guide's
      # denominators and its pseudo-gene grouping come from, and it matches the
      # guide ids examples/schmidt_tcell_barcs.R writes.
      sgRNA = paste0(set$guide, "|", set_tag),
      Gene = set$gene,
      normalised,
      check.names = FALSE
    )
  })
  do.call(rbind, blocks)
}

# Run `mageck test` on a prepared count table with the published options.
# `mageck test` shells out to the RRA binary by bare name, so an out-of-PATH
# install (a conda environment, say) needs its bin directory on PATH too.
schmidt_run_mageck <- function(table, controls, screen_dir, prefix,
                               mageck = Sys.getenv("MAGECK", unset = "mageck")) {
  if (nchar(Sys.which(mageck)) == 0L && !file.exists(mageck)) {
    stop("Cannot find the `mageck` executable. Install MAGeCK and set MAGECK=<path>.",
         call. = FALSE)
  }
  if (file.exists(mageck)) {
    Sys.setenv(PATH = paste(dirname(normalizePath(mageck)), Sys.getenv("PATH"), sep = ":"))
  }
  dir.create(screen_dir, recursive = TRUE, showWarnings = FALSE)

  counts_path <- file.path(screen_dir, paste0(prefix, "_counts.txt"))
  controls_path <- file.path(screen_dir, paste0(prefix, "_control_sgrna.txt"))
  log_path <- file.path(screen_dir, paste0(prefix, ".mageck.log"))
  write.table(table, counts_path, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(controls, controls_path)

  # --paired matches -t and -c positionally, so both lists are in donor order.
  donors <- sort(unique(sub(
    "_(low|high)$", "", setdiff(colnames(table), c("sgRNA", "Gene"))
  )))
  status <- system2(
    mageck,
    c("test", "-k", shQuote(counts_path),
      "-t", paste0(donors, "_high", collapse = ","),
      "-c", paste0(donors, "_low", collapse = ","),
      "--paired", "--norm-method", "none",
      "--control-sgrna", shQuote(controls_path),
      "-n", shQuote(file.path(screen_dir, prefix))),
    stdout = log_path, stderr = log_path
  )
  if (!identical(status, 0L)) {
    stop("mageck test failed for ", prefix, "; see ", log_path, call. = FALSE)
  }
  read.delim(file.path(screen_dir, paste0(prefix, ".gene_summary.txt")),
             check.names = FALSE)
}

# The published hit rule, quoted from the Methods: "Gene hits were classified as
# having a median absolute log2-fold change >0.5 and a false discovery rate
# (FDR) <0.05". RRA is one-sided in each direction, so the single FDR per gene
# is the smaller of the two tails.
schmidt_mageck_hits <- function(genes) {
  abs(genes[["pos|lfc"]]) > 0.5 &
    pmin(genes[["pos|fdr"]], genes[["neg|fdr"]]) < 0.05
}
