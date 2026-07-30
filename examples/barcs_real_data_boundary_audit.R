#!/usr/bin/env Rscript

# Fraction of converged guide fits at the rho-hat = 0 boundary in the three
# real-data analyses used to discuss finite-sample calibration.

options(stringsAsFactors = FALSE)

inputs <- list(
  `HT-29` = "results/chronos_tzelepis/beta_binomial_guide_results.csv.gz",
  IL2RA = "results/waterbear_facs/beta_binomial_all_bins_guide_results.csv.gz",
  Liang = Sys.glob(
    "results/liang_cas13/*_longitudinal_barcs_guide.csv.gz"
  )
)

pieces <- lapply(names(inputs), function(dataset) {
  files <- inputs[[dataset]]
  if (!length(files) || any(!file.exists(files))) {
    stop("Missing BARCS guide results for ", dataset)
  }
  result <- do.call(rbind, lapply(files, function(path) {
    read.csv(gzfile(path))
  }))
  usable <- result$converged & is.finite(result$rho)
  data.frame(
    dataset = dataset,
    guides = nrow(result),
    converged_guides = sum(usable),
    rho_zero_guides = sum(result$rho[usable] == 0),
    rho_zero_fraction = mean(result$rho[usable] == 0),
    median_rho = median(result$rho[usable]),
    median_mean_cpm = median(result$mean_cpm[usable])
  )
})

audit <- do.call(rbind, pieces)
write.csv(
  audit,
  "data/derived/barcs_real_data_rho_boundary.csv",
  row.names = FALSE
)
print(audit, row.names = FALSE)
