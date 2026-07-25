# Liang HAP1 full real-data case

Source study: Liang et al. transcriptome-scale Cas13 fitness screens
(Cell Genomics, 2026; BioProject PRJNA1344834).

## Full processed-count analysis

`liang-hap1-full-counts.csv` contains all 56,174 HAP1 guides for which
Liang deposited processed endpoint values. The four columns are the two
day-0 and two day-14 libraries. Values are the deposited normalized,
ComBat-corrected values rounded to integer pseudo-counts, exactly as in
the manuscript sensitivity analysis. `liang-hap1-metadata.csv` retains
the corresponding complete-library totals. Non-targeting guides are
assigned to four-guide pseudo-genes only for gene-level null calibration;
their guide identities and counts are unchanged.

This is a complete real screen, not a selected-gene toy case. Because the
deposited values have already been normalized and corrected, it remains a
processed-count sensitivity analysis rather than a raw-count likelihood
analysis.

## Deposited FASTQ analysis

`liang-hap1-real-fastq-manifest.csv` lists the four real HAP1 endpoint
FASTQs from ENA, including run accessions, HTTPS URLs, MD5 checksums, and
compressed byte sizes. Together they are about 1.7 GB.
`liang-hap1-full-guide-library.csv` contains all 56,322 library guides.
For browser analysis, rename each downloaded file to the
`browser_filename` in the manifest so it matches the sample metadata.
From the repository root, run:

```sh
bash scripts/run_liang_hap1_real_case.sh
```

The runner downloads one deposited FASTQ at a time with resume and
published-MD5 verification, applies the article-matched Cutadapt/Bowtie
workflow, writes a BARCS-ready raw count matrix, and deletes the temporary
compressed read file after counting.

## Small validation fixture

`liang-hap1-counts.csv` contains 72 real Liang guide rows: four guides
each for four essential/protein-coding genes, two published lncRNA
signals, two TPM-zero lncRNA nulls, and 40 non-targeting controls. Values
are Liang's deposited normalized, ComBat-corrected processed counts rounded
to integer pseudo-counts. `liang-hap1-metadata.csv` preserves full-library
totals calculated before the 72-guide subset was selected.

This is a compact same-input sensitivity example. It is not a raw-count
likelihood analysis and should not be used to replace the complete benchmark.

The four small `.fastq.gz` files are synthetic software-test fixtures, not
a biological example and not Liang's
deposited sequencing reads. They use the real selected 23-nt Liang guide
sequences and contain one exact-matching read per rounded selected-guide
pseudo-count. Ten percent of reads contain the reverse-complement spacer to
exercise orientation detection. Alignment should reproduce
`liang-hap1-fastq-expected-counts.csv` exactly.

Use `liang-hap1-guide-library.csv` as the guide library and
`liang-hap1-fastq-metadata.csv` as sample metadata. Fit `day14` while
adjusting for `replicate`.

The FASTQ-derived analysis uses totals from the selected synthetic library;
therefore its estimates need not equal the count example, which retains the
complete Liang library denominators.
