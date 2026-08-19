# Liang et al. Cas13 benchmark sources

Inputs for `examples/liang_cas13_benchmark.R`, the transcriptome-scale
RfxCas13d fitness screens of Liang et al. (Cell Genomics, 2026; BioProject
[PRJNA1344834](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1344834)).

Nothing in this directory is versioned. One script fetches and prepares all of
it:

```sh
Rscript scripts/prepare_liang_cas13.R
```

It needs the `readxl` package, downloads only what is missing, and leaves a
prepared tree of about 110 MB.

## Downloaded article supplements

The two supplementary workbooks are static files attached to the article:

| Local file | Direct download |
|---|---|
| `liang_table_s1.xlsx` | `https://ars.els-cdn.com/content/image/1-s2.0-S2666979X26001151-mmc2.xlsx` |
| `liang_table_s2.xlsx` | `https://ars.els-cdn.com/content/image/1-s2.0-S2666979X26001151-mmc3.xlsx` |

SHA-256:

```text
df55475c7e052009352228c14a56387858c3f60d69965808c69dff8613aa95ec  liang_table_s1.xlsx
27305b20b957c309d11434fe77c72b3a949990ad918edb2d54f03b3c6d80aed6  liang_table_s2.xlsx
```

## Downloaded run metadata

Two live API queries supply the sequencing-run manifest. These are generated
per request and their bytes change as SRA and ENA revise submission metadata,
so they are deliberately not checksummed:

- `PRJNA1344834_runinfo.csv` —
  `https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA1344834`
- `PRJNA1344834_ena_fastq.tsv` —
  `https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1344834&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes&format=tsv`

The prepare script asserts the run manifest still resolves to exactly 30
day-0/day-7/day-14 amplicon runs with one ENA FASTQ URL each, so a metadata
change that breaks the join fails loudly rather than silently.

## Prepared by `scripts/prepare_liang_cas13.R`

| Local file | Extracted from |
|---|---|
| `guide_library.tsv` | Table S1 sheet `S1B` — 56,322 unique 23-nt guides, asserted |
| `guide_reference.fa` | The same library as FASTA, for `bowtie-build` |
| `lncrna_expression.tsv` | Table S1 sheets `S1C` (total RNA) and `S1E` (mRNA) |
| `published_liang_rra.tsv` | Table S2 sheets `S2F`–`S2J` — published day-7 and day-14 RRA calls |
| `published_processed_counts_HAP1.tsv` etc. | Table S2 sheets `S2A`–`S2E` — deposited normalized counts |
| `endpoint_sample_manifest.tsv` | Run info joined to ENA FASTQ URLs, 30 rows |

`lncrna_expression.tsv` holds the per-cell-line maximum of the total-RNA and
mRNA TPM. The study's screen filter is TPM > 0, and taking the maximum of the
two assays avoids calling a gene unexpressed only because one library type did
not detect it. The benchmark uses TPM == 0 lncRNAs as biological null controls.

Sheets cover five cell lines — HAP1, HEK293FT, K562, MDA-MB-231, and THP1 —
and the prepare script keeps all five. The benchmark analyzes four: K562 is
excluded there because its deposited day-0 samples do not line up with the
longitudinal design.

The deposited counts are fractional, having been normalized, ComBat-corrected,
and outlier-processed. `examples/liang_cas13_benchmark.R` rounds them once and
hands the identical pseudo-count matrix to every method, which makes the Liang
comparison a same-input processed-count sensitivity analysis rather than a
likelihood-faithful raw-count one. The rounding magnitude is recorded in
`data/derived/liang_cas13_input_audit.csv`.

## Optional: recounting from FASTQ

The benchmark does not need this. It reproduces the article's own read
counting from the raw runs, and writes `counts/<sample>.counts.tsv` plus a
per-sample log:

```sh
bash scripts/queue_liang_cas13_counts.sh [parallelism]
```

That queues one restartable
[`scripts/count_liang_cas13_run.sh`](../../../scripts/count_liang_cas13_run.sh)
task per manifest row. Each streams its FASTQ, trims the published anchors
with Cutadapt, aligns to the Bowtie guide index, and tabulates guides. No
FASTQ is retained on disk, and finished samples are skipped on rerun.

Required: `pueue`, `bowtie` and `bowtie-build`, and Cutadapt at
`.venv/bin/cutadapt` (override with `CUTADAPT`, `BOWTIE`, `BOWTIE_BUILD`).
Reads stream through `fastq-dump` when sra-tools is installed, otherwise over
HTTPS from ENA. The `guide_reference.*.ebwt` index files are built on first
use.

The article's prose says "up to three mismatches" while its reported parameter
is `-v 1`; the counting script uses `-v 1` and records the discrepancy in each
sample log.
