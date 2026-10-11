# Liu et al. genome-wide T cell IFNγ screen — BARCS re-analysis

BARCS re-analysis of the genome-wide arm of the in vivo tumour-infiltrating
T cell IFNγ screen from Liu Q, Chen PA, Urs E, *et al.*, "In vivo genome-wide
CRISPR screens of human T cells in solid tumours," *Nature* (2026),
doi:10.1038/s41586-026-10906-9. The source is Supplementary Table 2, the
screen behind the paper's main Fig. 4.

This is the genome-wide companion to the focused sub-library re-analysis in
`results/liu_tcell/`. It is a re-analysis of published median-normalised
counts, not an independent replication.

## What differs from the sub-library

The genome-wide screen (19,113 genes, ~77k guides) is reported with **two donor
replicates** (`_r0`, `_r1`), not per-mouse counts. BARCS therefore fits a
donor-adjusted model, `~ gate + donor`, with the sort gate (IFNγ-high vs
IFNγ-low) as the tested term. There is no per-mouse structure to pair on, so the
`~ gate + mouse` sensitivity fit that is central to the sub-library analysis is
not available here.

## Layout

These tables are generated output and are **not** versioned — the genome-wide
count matrix (~4 MB) and per-guide table (~20 MB) exceed the compact-summary
policy in `results/README.md`. Regenerate them from the scripts below; the run
is about a minute.

- `input/` — `counts.csv` (guide × four libraries: two gates × two donors) and
  `metadata.csv` (gate, donor, library totals).
- `output/` — gene-level (`gate-donor_GENES.csv`) and guide-level
  (`gate-donor_GUIDES.csv`) BARCS results.

## Regenerate

From the repository root:

```
Rscript scripts/prepare_liu_genomewide.R    # Supp. Table 2 -> input/
Rscript examples/liu_genomewide_barcs.R     # input/ -> output/  (~1 min, multicore)
Rscript examples/liu_genomewide_publication_comparison.R  # figure vs MAGeCK
```

`prepare_liu_genomewide.R` reshapes the sgRNA-summary sheet (one row per
guide-and-donor) into the four-library matrix and rounds the normalised counts
half-up to integers; `data/raw/liu_tcell/README.md` records the archive
provenance. `liu_genomewide_barcs.R` fits `bb_screen(term = "gate")`, calibrates
against the 1,000 non-targeting (NTCTRL) guides with `bb_calibrate_controls()`,
and summarises genes with `bb_gene_consistency()`. The comparison script writes
`figures/liu_genomewide_publication_comparison.png` and `.pdf`.

## Notes

- Gene-level effects agree with the published MAGeCK log2 fold changes at
  Spearman ρ ≈ 0.82 over 19,105 shared genes, and IFNG is recovered as a strong
  depleted control.
- The guide-level non-targeting-control calibration is ×1.00 (no correction),
  as in the sub-library: the 1,000 NTCTRL guides give a well-behaved null
  (calibrated `p < 0.05` rate 0.003). The gene-level `null_scale` reported in the
  output (≈3.09) is the global median-absolute-deviation spread of all ~19k gene
  statistics, not a control-derived overdispersion correction — with only one
  control "gene" (NTCTRL), `bb_gene_consistency()` falls back to that global
  scale. It reflects the natural width of the genome-wide gene-statistic
  distribution, not a calibration cost of the two-donor design.
- The non-targeting controls sit at higher read abundance than targeting guides
  (median count ratio ≈2.44× here, versus ≈1.12× in the sub-library) — expected,
  since non-cutting guides do not drop out. This abundance inflation does **not**
  drive the calibration above: the guide-level control scale is 1.00 whether or
  not the controls are abundance-matched to the targeting guides, and the
  gene-level scale is global. It is noted as a property of the data, not a
  correction applied.
- BARCS effects are logit-scale beta-binomial coefficients, comparable to MAGeCK
  in sign and rank only, not magnitude. The `NTCTRL` control pseudo-gene pools
  1,000 guides, so its gene-level row has an anomalously small standard error
  and should not be read as a gene call; the guide-level control null is the
  relevant calibration check.
