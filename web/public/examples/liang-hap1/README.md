# Liang HAP1 BARCS examples

Source study: Liang et al. transcriptome-scale Cas13 fitness screens 
(Cell Genomics, 2026; BioProject PRJNA1344834).

## Count-matrix example

`liang-hap1-counts.csv` contains 72 real Liang guide rows: four guides 
each for four essential/protein-coding genes, two published lncRNA 
signals, two TPM-zero lncRNA nulls, and 40 non-targeting controls. Values 
are Liang's deposited normalized, ComBat-corrected processed counts rounded 
to integer pseudo-counts. `liang-hap1-metadata.csv` preserves full-library 
totals calculated before the 72-guide subset was selected.

This is a compact same-input sensitivity example. It is not a raw-count 
likelihood analysis and should not be used to replace the complete benchmark.

## FASTQ example

The four `.fastq.gz` files are synthetic teaching fixtures, not Liang's 
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
