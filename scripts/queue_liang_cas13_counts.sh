#!/usr/bin/env bash

# Submit restartable, streamed read-counting tasks to a dedicated pueue group.
# Usage: bash scripts/queue_liang_cas13_counts.sh [parallelism]

set -euo pipefail

parallelism=${1:-2}
project_root=$(cd "$(dirname "$0")/.." && pwd)
manifest="${project_root}/data/raw/liang_cas13/endpoint_sample_manifest.tsv"
output_dir="${project_root}/data/raw/liang_cas13/counts"
group="liang-cas13"

if [[ ! -s "${manifest}" ]]; then
  echo "Run Rscript scripts/prepare_liang_cas13.R first." >&2
  exit 1
fi
if ! command -v pueue >/dev/null 2>&1; then
  echo "pueue is required for this long-running batch." >&2
  exit 1
fi

mkdir -p "${output_dir}"
if ! pueue group --json | grep -q "\"${group}\""; then
  pueue group add "${group}"
fi
pueue parallel "${parallelism}" --group "${group}"

tail -n +2 "${manifest}" |
  while IFS=$'\t' read -r run sample spots size_mb sra_path \
    cell_line day replicate fastq_url; do
    count_path="${output_dir}/${sample}.counts.tsv"
    if [[ -s "${count_path}" ]]; then
      echo "${sample}: already counted"
      continue
    fi
    pueue add \
      --group "${group}" \
      --label "${sample}" \
      -w "${project_root}" \
      -- bash scripts/count_liang_cas13_run.sh \
        "${sample}" "${run}" "${fastq_url}" \
        "data/raw/liang_cas13/counts"
  done

pueue status --group "${group}"
