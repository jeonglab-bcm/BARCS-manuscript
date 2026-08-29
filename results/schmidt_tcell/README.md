# Schmidt et al. T cell CRISPRa/CRISPRi screens — BARCS versus MAGeCK

A controlled head-to-head on the genome-wide CRISPRa and CRISPRi screens of
Schmidt R, Steinhart Z, Layeghi M, *et al.*, "CRISPR activation and interference
screens decode stimulation responses in primary human T cells," *Science* 376,
eabj4008 (2022), doi:10.1126/science.abj4008, from the raw sgRNA read counts in
GEO GSE174255.

Both methods are run here on the identical raw counts — BARCS, and the paper's
own MAGeCK recipe — so every difference is the model rather than the input.
`BARCS_summary.md` is the full write-up: fit diagnostics, the reproduced
published biology, the head-to-head, and the held-out null calibration.

**The short version.** The two methods order the screens much the same way
(Spearman 0.84–0.85) and recover the paper's biology equally well: BARCS leads
inside the top 50 reference genes (36 vs 34 of 43), MAGeCK by one gene from the
top 100 outward. They differ in calibration — on 660 held-out nontargeting
pseudo-genes BARCS runs at 0.52× its nominal null rate and MAGeCK at 1.24× of the
matched expectation, so BARCS's shorter hit lists are conservatism rather than
weakness, with power still recoverable.

## Layout

- `comparison/` — the committed tables. `method_concordance.csv` (per-screen
  agreement, null scales and hit counts under three rules), `null_calibration.csv`
  (the held-out check, one row per screen × fold × method),
  `positive_control_panel.csv` (rank and call of each reference gene under each
  method), `tcr_pathway_precision.csv` (KEGG TCR-signalling genes among each
  method's top N), and `<screen>_hits.csv` (every gene either method calls, plus
  the reference panel, with both methods' numbers side by side).
- `input/`, `output/`, `mageck/`, `mageck_holdout/` — regenerated, not versioned.

## Regenerate

From the repository root:

```
Rscript scripts/prepare_schmidt_tcell.R              # GEO GSE174255 -> input/
Rscript examples/schmidt_tcell_barcs.R               # input/ -> output/
Rscript examples/schmidt_tcell_mageck.R              # input/ -> mageck/
Rscript examples/schmidt_tcell_null_calibration.R    # -> null_calibration.csv
Rscript examples/schmidt_tcell_method_comparison.R   # -> comparison/, figure
```

Steps 3 and 4 need the `mageck` executable; set `MAGECK` to its path if it is not
on `PATH`. `BARCS_NCORES` sets the guide-fitting parallelism (default 4); the
eight library-set fits take about 70 seconds at 8 cores. `R/schmidt_tcell.R`
holds the pieces the scripts share — the library sets, the nontargeting
pseudo-gene construction, the held-out fold split, and the MAGeCK table and call.

## The models

BARCS fits one guide-level beta-binomial regression per library set over all
twelve sorted and unsorted libraries,

```
~ donor + assay * bin        bin in {low, unsorted, high}, assay in {IFNG, IL2}
```

then moderates the per-guide dispersions across the guides sharing that design
(5 residual df → ~140 in CRISPRa, ~107 in CRISPRi), calibrates on the
nontargeting controls, and summarises genes. Both cytokine contrasts come off the
single fit. The two library sets are never rescaled onto each other; they meet
only at the gene level, where they supply 3 + 3 = 6 guides.

The published pipeline normalises every library to its total read count, merges
the two sets into one table, and runs `mageck test --norm-method none --paired
--control-sgrna`, pairing by donor.

Dispersion moderation is applied here on the strength of the held-out check, not
as a general default: on the 268-guide Liu screen in `results/liu_tcell/` the
same step costs four published hits and raises the null scale.

## Figure

`examples/schmidt_tcell_method_comparison.R` writes
`figures/schmidt_tcell_method_comparison.{png,pdf}` (figures are generated output
and are not versioned): effect concordance per screen, where each method ranks
the genes the paper names, hits per screen under three rules, the held-out null
calibration, and the hit-count spread across held-out fits.
