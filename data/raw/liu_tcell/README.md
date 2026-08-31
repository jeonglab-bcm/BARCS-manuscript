# Liu et al. in vivo T cell screen sources

Source for the BARCS re-analysis in `results/liu_tcell/`. Liu Q, Chen PA, Urs E,
*et al.*, "In vivo genome-wide CRISPR screens of human T cells in solid
tumours," *Nature* (2026),
[doi:10.1038/s41586-026-10906-9](https://doi.org/10.1038/s41586-026-10906-9).

Nothing in this directory is versioned. One script fetches and prepares it:

```sh
Rscript scripts/prepare_liu_tcell.R
```

It needs the `readxl` package, downloads only what is missing, and writes the
four model inputs into `results/liu_tcell/input/`.

## Downloaded article supplement

The supplementary workbooks are a single archive attached to the article:

- `41586_2026_10906_MOESM3_ESM.zip` — supplementary tables archive
  (`https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-026-10906-9/MediaObjects/41586_2026_10906_MOESM3_ESM.zip`),
  from which `Supplementary_Table_5.xlsx` (focused sub-library screen) and
  `Supplementary_Table_2.xlsx` (genome-wide screen, main Fig. 4) are extracted.

Two re-analyses draw from this archive: the focused sub-library
(`results/liu_tcell/`, from Supplementary Table 5, below) and the genome-wide
screen (`results/liu_genomewide/`, from Supplementary Table 2, via
`scripts/prepare_liu_genomewide.R`).

## What the script derives

Supplementary Table 5 is the only screen in the paper reporting per-mouse
counts. The re-analysis uses the two pooled sgRNA-summary sheets:

- Arm A (CD3-scFv vs A375low): `CD3scFv_A375low_3Donor_sgRNASum`
- Arm B (NY-ESO-1 TCR vs WT A375): `WT_A375_2donors_sgRNASum`

Each guide's `control_count` (IFNγ-low gate) and `treatment_count` (IFNγ-high
gate) hold the MAGeCK median-normalised counts as `/`-separated per-library
values. The script splits them, rounds half-up to integers, interleaves the two
gates per mouse, and derives the sample metadata (gate, donor, mouse, and the
per-library totals) from the column structure. The per-donor and per-mouse
sheets are separately renormalised and are not used.

The screen FASTQs were never deposited; only the study's RNA-seq is public
(GSE330227). These normalised counts are therefore the most raw form of the
screen obtainable, so the analysis is a re-analysis of published normalised
counts, not an independent replication.
