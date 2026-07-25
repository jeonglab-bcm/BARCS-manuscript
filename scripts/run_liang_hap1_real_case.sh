#!/usr/bin/env bash

# Reproduce the complete Liang et al. HAP1 day-0/day-14 raw-read input.
#
# The four deposited ENA FASTQs (~1.7 GB compressed) are streamed through the
# article-matched Cutadapt/Bowtie workflow. FASTQs are not retained. Restarting
# is safe because completed sample count files are skipped.

set -euo pipefail

project_root=$(cd "$(dirname "$0")/.." && pwd)
raw_dir="${project_root}/data/raw/liang_cas13"
manifest="${raw_dir}/endpoint_sample_manifest.tsv"
count_dir="${raw_dir}/counts"
output_dir="${project_root}/results/liang_hap1_raw"

if [[ ! -s "${manifest}" || ! -s "${raw_dir}/guide_library.tsv" ]]; then
  echo "Run Rscript scripts/prepare_liang_cas13.R first." >&2
  exit 1
fi

mkdir -p "${count_dir}" "${output_dir}"
cd "${project_root}"

while IFS=$'\t' read -r run sample spots size_mb sra_path \
  cell_line day replicate fastq_url; do
  if [[ "${cell_line}" != "HAP1" ]]; then
    continue
  fi
  if [[ "${day}" != "0" && "${day}" != "14" ]]; then
    continue
  fi
  bash scripts/count_liang_cas13_run.sh \
    "${sample}" "${run}" "${fastq_url}" \
    "data/raw/liang_cas13/counts"
done < <(tail -n +2 "${manifest}")

Rscript - "${count_dir}" "${output_dir}" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
count_dir <- args[[1L]]
output_dir <- args[[2L]]
guide <- read.delim(
  file.path("data", "raw", "liang_cas13", "guide_library.tsv"),
  check.names = FALSE
)
samples <- c(
  "HAP1_Day00_R1", "HAP1_Day00_R2",
  "HAP1_Day14_R1", "HAP1_Day14_R2"
)

counts <- matrix(
  0L,
  nrow = nrow(guide),
  ncol = length(samples),
  dimnames = list(guide$sgrna, samples)
)
for (column in seq_along(samples)) {
  path <- file.path(count_dir, paste0(samples[[column]], ".counts.tsv"))
  if (!file.exists(path)) {
    stop("Missing completed count file: ", path)
  }
  observed <- read.delim(path, check.names = FALSE)
  index <- match(observed$sgrna, guide$sgrna)
  if (anyNA(index)) {
    stop("Aligned guide absent from Liang library in ", path)
  }
  counts[index, column] <- as.integer(observed$count)
}

control <- guide$target_group == "non-targeting"
gene <- guide$gene
control_index <- which(control)
gene[control_index] <- sprintf(
  "NTC_%03d",
  ceiling(seq_along(control_index) / 4)
)
count_table <- data.frame(
  guide = guide$sgrna,
  gene = gene,
  control = tolower(as.character(control)),
  counts,
  check.names = FALSE
)
metadata <- data.frame(
  sample = samples,
  day14 = c(0L, 0L, 1L, 1L),
  replicate = c("R1", "R2", "R1", "R2"),
  total = as.integer(colSums(counts))
)
write.csv(
  count_table,
  file.path(output_dir, "liang-hap1-raw-counts.csv"),
  row.names = FALSE,
  quote = FALSE
)
write.csv(
  metadata,
  file.path(output_dir, "liang-hap1-raw-metadata.csv"),
  row.names = FALSE,
  quote = FALSE
)
message(
  "Wrote ", nrow(count_table), " guides across ", ncol(counts),
  " real libraries to ", output_dir, "."
)
RS

echo
echo "Real Liang HAP1 inputs are ready:"
echo "  ${output_dir}/liang-hap1-raw-counts.csv"
echo "  ${output_dir}/liang-hap1-raw-metadata.csv"
