#!/usr/bin/env Rscript

# External validation of method-specific GSE70038 top-200 hits against the
# Hart reference-essential set used in the independent A375 benchmark.

options(stringsAsFactors = FALSE)

load(file.path("CB2", "data", "Sanson_CRISPRn_A375.rda"))
reference_essential <- unique(Sanson_CRISPRn_A375$egenes)
barcs <- read.csv("results/gse70038/beta_binomial_gene_results.csv")
mageck <- read.csv("results/gse70038/mageck_gene_results.csv")
terminal_terms <- c(
  "GSC0131_end", "GSC0827_end", "NSCCB660_end", "NSCU5_end"
)
testable_universe <- intersect(unique(barcs$gene), unique(mageck$gene))

top_depleted <- function(result, terminal_term, number = 200L) {
  x <- result[
    result$term == terminal_term &
      is.finite(result$estimate) &
      is.finite(result$p_value),
  ]
  x <- x[order(x$estimate, x$p_value, x$gene), ]
  unique(head(x$gene, number))
}

pieces <- lapply(terminal_terms, function(terminal_term) {
  barcs_top <- top_depleted(barcs, terminal_term)
  mageck_top <- top_depleted(mageck, terminal_term)
  sets <- list(
    BARCS = setdiff(barcs_top, mageck_top),
    `MAGeCK-MLE` = setdiff(mageck_top, barcs_top)
  )
  do.call(rbind, lapply(names(sets), function(method) {
    selected <- sets[[method]]
    not_selected <- setdiff(testable_universe, selected)
    contingency <- matrix(
      c(
        sum(selected %in% reference_essential),
        sum(!selected %in% reference_essential),
        sum(not_selected %in% reference_essential),
        sum(!not_selected %in% reference_essential)
      ),
      nrow = 2,
      byrow = TRUE
    )
    test <- fisher.test(contingency, alternative = "greater")
    data.frame(
      coefficient = sub("_end$", "", terminal_term),
      method = method,
      exclusive_genes = length(selected),
      reference_essential = sum(selected %in% reference_essential),
      reference_essential_fraction = mean(selected %in% reference_essential),
      enrichment_odds_ratio = unname(test$estimate),
      one_sided_p_value = test$p.value
    )
  }))
})

result <- do.call(rbind, pieces)
write.csv(
  result,
  "data/derived/gse70038_reference_essential_enrichment.csv",
  row.names = FALSE
)
print(result, row.names = FALSE)
