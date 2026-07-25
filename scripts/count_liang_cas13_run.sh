#!/usr/bin/env bash

# Stream one Liang et al. SRA/ENA FASTQ through the published anchor-trimming
# and Bowtie guide-alignment workflow.  No FASTQ is retained on disk.

set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "usage: $0 SAMPLE RUN FASTQ_URL OUTPUT_DIR" >&2
  exit 2
fi

sample=$1
run=$2
fastq_url=$3
output_dir=$4
project_root=$(cd "$(dirname "$0")/.." && pwd)
reference_dir="${project_root}/data/raw/liang_cas13"
reference_fasta="${reference_dir}/guide_reference.fa"
index_prefix="${reference_dir}/guide_reference"
count_path="${output_dir}/${sample}.counts.tsv"
log_path="${output_dir}/${sample}.log"
temporary_path="${count_path}.partial"

cutadapt_executable=${CUTADAPT:-"${project_root}/.venv/bin/cutadapt"}
bowtie_executable=${BOWTIE:-bowtie}
bowtie_build_executable=${BOWTIE_BUILD:-bowtie-build}

for executable in "${cutadapt_executable}" \
  "${bowtie_executable}" "${bowtie_build_executable}"; do
  if ! command -v "${executable}" >/dev/null 2>&1; then
    echo "Required executable is missing: ${executable}" >&2
    exit 1
  fi
done
if [[ ! -s "${reference_fasta}" ]]; then
  echo "Run scripts/prepare_liang_cas13.R before counting reads." >&2
  exit 1
fi

mkdir -p "${output_dir}"
if [[ -s "${count_path}" ]]; then
  echo "${sample}: count file already exists; skipping."
  exit 0
fi

if [[ ! -s "${index_prefix}.1.ebwt" ]]; then
  "${bowtie_build_executable}" "${reference_fasta}" "${index_prefix}"
fi

trap 'rm -f "${temporary_path}"' EXIT
{
  echo "sample=${sample}"
  echo "run=${run}"
  echo "url=${fastq_url}"
  echo "started_utc=$(date -u +%FT%TZ)"
  echo "cutadapt_version=$("${cutadapt_executable}" --version)"
  echo "bowtie_version=$("${bowtie_executable}" --version | head -1)"
  if command -v fastq-dump >/dev/null 2>&1; then
    echo "transport=fastq-dump remote stream"
  else
    echo "transport=ENA HTTPS gzip stream"
  fi
} >"${log_path}"

stream_fastq() {
  if command -v fastq-dump >/dev/null 2>&1; then
    fastq-dump --stdout "${run}"
  else
    for executable in curl gzip; do
      if ! command -v "${executable}" >/dev/null 2>&1; then
        echo "Fallback transport requires ${executable}." >&2
        return 1
      fi
    done
    curl --location --fail --retry 5 --retry-all-errors --silent \
      --show-error "${fastq_url}" |
      gzip -dc
  fi
}

# The two Cutadapt invocations reproduce the article's stated anchor rules.
# The paper says "up to three mismatches" in prose but reports `-v 1`; we use
# the executable parameter (`-v 1`) and preserve this discrepancy in the log.
stream_fastq |
  "${cutadapt_executable}" \
    -g CTGGTCGGGGTTTGAAAC -e 0.2 -O 5 --discard-untrimmed - 2>>"${log_path}" |
  "${cutadapt_executable}" \
    -a TTTTTGAATTCGCTAGCT -e 0.1 -O 5 --minimum-length 15 \
    --discard-untrimmed - 2>>"${log_path}" |
  "${bowtie_executable}" -v 1 -m 3 --best -q \
    "${index_prefix}" - --suppress 1,2,4,5,6,7,8 2>>"${log_path}" |
  LC_ALL=C sort |
  uniq -c |
  awk 'BEGIN {OFS="\\t"; print "sgrna", "count"} {print $2, $1}' \
    >"${temporary_path}"

mv "${temporary_path}" "${count_path}"
echo "finished_utc=$(date -u +%FT%TZ)" >>"${log_path}"
echo "${sample}: wrote ${count_path}"
