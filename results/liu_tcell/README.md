# Liu et al. in vivo T cell screen — BARCS re-analysis

BARCS re-analysis of the focused sub-library screen from Liu Q, Chen PA, Urs E,
*et al.*, "In vivo genome-wide CRISPR screens of human T cells in solid
tumours," *Nature* (2026), doi:10.1038/s41586-026-10906-9 (Supplementary
Table 5 — the only screen in the paper publishing per-mouse counts).

`BARCS_PI_summary.md` is the full write-up: fit diagnostics, the reproduced
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

From the repository root:

```
Rscript examples/liu_tcell_barcs.R
```

The script reads `input/` and rewrites `output/` with the BARCS R package
directly (no BARCS Studio). For each arm it fits `bb_screen(term = "gate")`,
calibrates against the non-targeting controls with `bb_calibrate_controls()`,
and summarises genes with `bb_gene_consistency()`.
