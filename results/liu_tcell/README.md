# Liu et al. in vivo T cell screen — BARCS re-analysis

BARCS re-analysis of the focused sub-library screen from Liu Q, Chen PA, Urs E,
*et al.*, "In vivo genome-wide CRISPR screens of human T cells in solid
tumours," *Nature* (2026), doi:10.1038/s41586-026-10906-9 (Supplementary
Table 5 — the only screen in the paper publishing per-mouse counts).

`BARCS_summary.md` is the full write-up: fit diagnostics, the reproduced
published hits, the donor-versus-mouse pairing result, and the limits. This is
a re-analysis of published normalised counts, not an independent replication.

## Layout

- `input/` — model inputs, one guide-by-library count matrix and one metadata
  table per arm (`armA_*`, `armB_*`). Counts are MAGeCK median-normalised,
  integer-valued, and each library's `total` equals its column sum.
- `output/` — BARCS results, gene-level (`*_GENES.csv`, 13 cols) and
  guide-level (`*_GUIDES.csv`, 17 cols), for each arm and pairing factor:
  `armA_gate-donor`, `armA_gate-mouse` (Result 3 sensitivity fit), and
  `armB_gate-donor`.

## Regenerate

From the repository root, two steps run the analysis from the published archive:

```
Rscript scripts/prepare_liu_tcell.R   # download Supp. Table 5 -> input/
Rscript examples/liu_tcell_barcs.R    # input/ -> output/
```

`prepare_liu_tcell.R` downloads the Nature supplementary archive, parses
Supplementary Table 5, and writes the four files in `input/` (see
`data/raw/liu_tcell/README.md` for the sheet-level provenance). The inputs are
also committed, so the second step alone regenerates the outputs on a clean
checkout. `liu_tcell_barcs.R` reads `input/` and rewrites `output/` with the
BARCS R package directly (no BARCS Studio): for each arm it fits
`bb_screen(term = "gate")`, calibrates against the non-targeting controls with
`bb_calibrate_controls()`, and summarises genes with `bb_gene_consistency()`.

## Comparison with the publication

`examples/liu_tcell_publication_comparison.R` draws BARCS against the published
MAGeCK analysis (Extended Data Fig. 5) and writes
`figures/liu_tcell_publication_comparison.png` (figures are generated output and
are not versioned). It needs `readxl` and the prepared Supplementary Table 5:

```
Rscript examples/liu_tcell_publication_comparison.R
```

Panels: (a, b) BARCS effect versus the published MAGeCK gene log2 fold change
for each arm (Spearman ρ ≈ 0.90 / 0.92); (c) BARCS Arm A versus Arm B, the
cross-model agreement the paper shows with RRA scores in panel g; (d) a BARCS
Arm A volcano, the counterpart of the MAGeCK RRA volcano in panel f. BARCS
reports a logit-scale coefficient and MAGeCK a log2 fold change, so the panels
are read as concordance in sign and rank, not as a fitted slope.
