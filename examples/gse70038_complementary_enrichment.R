#!/usr/bin/env Rscript

# Functional follow-up of method-specific GSE70038 depletion hits.
#
# For each terminal coefficient, this script selects the 200 most depleted
# genes from BARCS and MAGeCK-MLE, partitions them into shared and
# method-specific sets, and submits the method-specific sets to Enrichr.
# The fixed top-200 rule gives the two methods equally sized inputs and avoids
# interpreting their differently calibrated FDR thresholds as comparable.

required_packages <- c("httr", "jsonlite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the required R packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

result_dir <- file.path("results", "gse70038")
derived_dir <- file.path("data", "derived")
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)

barcs <- read.csv(
  file.path(result_dir, "beta_binomial_gene_results.csv"),
  stringsAsFactors = FALSE
)
mageck <- read.csv(
  file.path(result_dir, "mageck_gene_results.csv"),
  stringsAsFactors = FALSE
)

terminal_terms <- c(
  "GSC0131_end",
  "GSC0827_end",
  "NSCCB660_end",
  "NSCU5_end"
)
method_tables <- list(BARCS = barcs, `MAGeCK-MLE` = mageck)

top_depleted <- function(result, terminal_term, number = 200L) {
  rows <- result[
    result$term == terminal_term &
      is.finite(result$estimate) &
      is.finite(result$p_value),
    ,
    drop = FALSE
  ]
  rows <- rows[order(rows$estimate, rows$p_value, rows$gene), ]
  unique(head(rows$gene, number))
}

partition_one_term <- function(terminal_term) {
  barcs_top <- top_depleted(method_tables$BARCS, terminal_term)
  mageck_top <- top_depleted(method_tables$`MAGeCK-MLE`, terminal_term)
  shared <- intersect(barcs_top, mageck_top)

  list(
    counts = data.frame(
      terminal_coefficient = sub("_end$", "", terminal_term),
      barcs_top = length(barcs_top),
      mageck_top = length(mageck_top),
      shared = length(shared),
      barcs_only = length(setdiff(barcs_top, mageck_top)),
      mageck_only = length(setdiff(mageck_top, barcs_top)),
      jaccard = length(shared) / length(union(barcs_top, mageck_top)),
      stringsAsFactors = FALSE
    ),
    lists = list(
      BARCS = setdiff(barcs_top, mageck_top),
      `MAGeCK-MLE` = setdiff(mageck_top, barcs_top)
    )
  )
}

partitions <- lapply(terminal_terms, partition_one_term)
names(partitions) <- terminal_terms
count_summary <- do.call(rbind, lapply(partitions, `[[`, "counts"))
write.csv(
  count_summary,
  file.path(derived_dir, "gse70038_complementary_hit_counts.csv"),
  row.names = FALSE
)

enrichr_base <- "https://maayanlab.cloud/Enrichr"
gene_set_libraries <- c(
  "GO_Biological_Process_2023",
  "Reactome_2022"
)

# Keep temporary input-file handling inside this wrapper so it is always
# removed after the request.
submit_gene_list <- function(genes, description) {
  gene_file <- tempfile(fileext = ".txt")
  on.exit(unlink(gene_file), add = TRUE)
  writeLines(genes, gene_file)
  response <- httr::RETRY(
    "POST",
    paste0(enrichr_base, "/addList"),
    body = list(
      list = httr::upload_file(gene_file),
      description = description
    ),
    encode = "multipart",
    times = 6L,
    pause_base = 2,
    pause_cap = 20,
    quiet = FALSE
  )
  httr::stop_for_status(response)
  httr::content(response, as = "parsed", type = "application/json")$userListId
}

read_enrichment <- function(user_list_id, library_name) {
  response <- httr::RETRY(
    "GET",
    paste0(enrichr_base, "/enrich"),
    query = list(
      userListId = user_list_id,
      backgroundType = library_name
    ),
    times = 6L,
    pause_base = 2,
    pause_cap = 20,
    quiet = FALSE
  )
  httr::stop_for_status(response)
  payload <- httr::content(response, as = "text", encoding = "UTF-8")
  records <- jsonlite::fromJSON(payload, simplifyVector = FALSE)[[library_name]]
  if (!length(records)) {
    return(data.frame())
  }

  do.call(rbind, lapply(records, function(record) {
    data.frame(
      rank = as.integer(record[[1L]]),
      term = as.character(record[[2L]]),
      p_value = as.numeric(record[[3L]]),
      odds_ratio = as.numeric(record[[4L]]),
      combined_score = as.numeric(record[[5L]]),
      overlapping_genes = paste(unlist(record[[6L]]), collapse = ";"),
      adjusted_p_value = as.numeric(record[[7L]]),
      stringsAsFactors = FALSE
    )
  }))
}

enrichment_results <- list()
result_index <- 0L
for (terminal_term in terminal_terms) {
  for (method_name in names(partitions[[terminal_term]]$lists)) {
    genes <- partitions[[terminal_term]]$lists[[method_name]]
    description <- paste(
      "GSE70038",
      sub("_end$", "", terminal_term),
      method_name,
      "top-200-specific depleted genes"
    )
    user_list_id <- submit_gene_list(genes, description)

    for (library_name in gene_set_libraries) {
      result <- read_enrichment(user_list_id, library_name)
      if (!nrow(result)) {
        next
      }
      result_index <- result_index + 1L
      result$terminal_coefficient <- sub("_end$", "", terminal_term)
      result$method_specific_set <- method_name
      result$gene_set_library <- library_name
      result$input_genes <- length(genes)
      enrichment_results[[result_index]] <- result
    }
    Sys.sleep(1)
  }
}

enrichment_table <- do.call(rbind, enrichment_results)
enrichment_table <- enrichment_table[
  order(
    enrichment_table$terminal_coefficient,
    enrichment_table$method_specific_set,
    enrichment_table$gene_set_library,
    enrichment_table$adjusted_p_value,
    enrichment_table$rank
  ),
]
write.csv(
  enrichment_table,
  file.path(derived_dir, "gse70038_complementary_enrichr.csv"),
  row.names = FALSE
)

significant_summary <- enrichment_table[
  enrichment_table$adjusted_p_value < 0.05,
]
group_id <- interaction(
  significant_summary$terminal_coefficient,
  significant_summary$method_specific_set,
  significant_summary$gene_set_library,
  drop = TRUE
)
significant_summary <- do.call(rbind, lapply(
  split(significant_summary, group_id),
  function(rows) head(rows, 3L)
))
write.csv(
  significant_summary,
  file.path(
    derived_dir,
    "gse70038_complementary_enrichr_significant_summary.csv"
  ),
  row.names = FALSE
)

print(count_summary, row.names = FALSE)
cat(
  "\nSaved ", nrow(enrichment_table), " enrichment records and ",
  nrow(significant_summary), " significant summary rows.\n",
  sep = ""
)
