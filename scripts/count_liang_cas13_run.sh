#!/usr/bin/env bash

# Download one Liang et al. ENA FASTQ with resume support, verify its published
# MD5, run the article-matched anchor-trimming/Bowtie workflow, then delete the
# temporary compressed file.

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
transport=${LIANG_FASTQ_TRANSPORT:-ena}
temporary_fastq=""

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

trap 'rm -f "${temporary_path}" "${temporary_fastq}"' EXIT
{
  echo "sample=${sample}"
  echo "run=${run}"
  echo "url=${fastq_url}"
  echo "started_utc=$(date -u +%FT%TZ)"
  echo "cutadapt_version=$("${cutadapt_executable}" --version)"
  echo "bowtie_version=$("${bowtie_executable}" --version | head -1)"
  if [[ "${transport}" == "sra" ]] &&
      command -v fastq-dump >/dev/null 2>&1; then
    echo "transport=fastq-dump remote stream"
  else
    echo "transport=ENA HTTPS resumable download with MD5 verification"
  fi
} >"${log_path}"

prepare_ena_fastq() {
  local ena_manifest="${reference_dir}/PRJNA1344834_ena_fastq.tsv"
  local expected_md5
  local actual_md5
  local complete=0
  if [[ ! -s "${ena_manifest}" ]]; then
    echo "Missing ENA checksum manifest: ${ena_manifest}" >&2
    return 1
  fi
  expected_md5=$(awk -F '\t' -v run="${run}" \
    'NR > 1 && $1 == run {print $3; exit}' "${ena_manifest}")
  if [[ ! "${expected_md5}" =~ ^[0-9a-fA-F]{32}$ ]]; then
    echo "No valid published MD5 was found for ${run}." >&2
    return 1
  fi
  temporary_fastq=$(mktemp "${TMPDIR:-/tmp}/barcs-${run}.fastq.gz.XXXXXX")
  if command -v aria2c >/dev/null 2>&1; then
    if aria2c \
        --allow-overwrite=true \
        --auto-file-renaming=false \
        --continue=true \
        --dir="$(dirname "${temporary_fastq}")" \
        --file-allocation=none \
        --max-connection-per-server=8 \
        --max-tries=10 \
        --min-split-size=5M \
        --out="$(basename "${temporary_fastq}")" \
        --retry-wait=2 \
        --split=8 \
        --summary-interval=30 \
        "${fastq_url}"; then
      complete=1
    fi
  else
    for attempt in 1 2 3 4 5; do
      if curl --location --fail --retry 5 --retry-all-errors \
          --continue-at - --output "${temporary_fastq}" \
          --silent --show-error "${fastq_url}"; then
        complete=1
        break
      fi
      echo "${sample}: ENA transfer attempt ${attempt} was interrupted; resuming." >&2
    done
  fi
  if [[ "${complete}" -ne 1 ]]; then
    echo "${sample}: ENA transfer did not complete after five resumptions." >&2
    return 1
  fi
  actual_md5=$(md5 -q "${temporary_fastq}")
  if [[ "${actual_md5}" != "${expected_md5}" ]]; then
    echo "${sample}: FASTQ MD5 mismatch (${actual_md5} != ${expected_md5})." >&2
    return 1
  fi
  echo "${sample}: verified ${run} (${actual_md5})." >&2
}

stream_fastq() {
  if [[ "${transport}" == "sra" ]] &&
      command -v fastq-dump >/dev/null 2>&1; then
    fastq-dump --stdout "${run}"
  else
    gzip -dc "${temporary_fastq}"
  fi
}

if [[ "${transport}" != "sra" ]] ||
    ! command -v fastq-dump >/dev/null 2>&1; then
  for executable in curl gzip md5; do
    if ! command -v "${executable}" >/dev/null 2>&1; then
      echo "ENA transport requires ${executable}." >&2
      exit 1
    fi
  done
  prepare_ena_fastq
fi

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
  awk 'BEGIN {OFS="\t"; print "sgrna", "count"} {print $2, $1}' \
    >"${temporary_path}"

if [[ "$(head -n 1 "${temporary_path}")" != $'sgrna\tcount' ]]; then
  echo "${sample}: generated count table is not valid tab-separated output." >&2
  exit 1
fi
mv "${temporary_path}" "${count_path}"
echo "finished_utc=$(date -u +%FT%TZ)" >>"${log_path}"
echo "${sample}: wrote ${count_path}"
