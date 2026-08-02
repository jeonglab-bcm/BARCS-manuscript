# BARCS: Beta-binomial Analysis and Regression for CRISPR Screens

CB² asks whether guide abundance differs between two groups. BARCS answers the
general coefficient question: does guide abundance follow dose, time, an
ordered phenotype, or an adjusted treatment effect?

BARCS keeps CB²'s sampling idea—binomial sequencing variation plus
between-library beta-binomial heterogeneity—but replaces the two-group mean
with an ordinary R model matrix. Guide counts remain the response; dose, time,
donor, batch, and interaction terms are predictors. Coefficients and named
contrasts use a Student t reference based on independent samples rather than
reads. BARCS is a standalone package; it is a generalisation of CB²'s
two-group test, not a replacement for that workflow.

## Install

The R package lives in
[`jeonglab-bcm/BARCS`](https://github.com/jeonglab-bcm/BARCS) and is tracked
here as the `BARCS/` submodule, so this repository pins the exact package
commit used for each manuscript and benchmark revision. CB2 is still tracked
as a second submodule, but only as a benchmark baseline and as the source of
the Sanson screen data — nothing in BARCS depends on it.

The manuscript LaTeX lives in the `overleaf/` submodule, which clones from the
public mirror [`jeonglab-bcm/BARCS-tex`](https://github.com/jeonglab-bcm/BARCS-tex)
so that the prose stays searchable and diffable on GitHub. Coauthors edit it in
Overleaf; `scripts/manuscript_sync.sh` keeps the mirror and Overleaf in step.

```sh
git clone --recurse-submodules https://github.com/jeonglab-bcm/BARCS-manuscript.git
cd BARCS-manuscript
```

Analysis scripts load the package through `R/load_barcs.R`, which prefers an
installed BARCS and otherwise loads the pinned submodule in place. A clean
checkout therefore needs no install step, though installing is faster:

```sh
R CMD INSTALL BARCS
```

## Use

`bbreg()` fits one guide, `bb_contrast()` tests a named contrast, and
`bb_screen()` runs a guide-by-guide screen with guide-to-gene aggregation.
A minimal longitudinal analysis with annotated counts, library totals, and
sample metadata:

```sh
Rscript examples/barcs_quickstart.R
```

## Layout

- `BARCS/` — pinned package submodule; the implementation and its own tests
- `R/load_barcs.R` — resolves the package for every analysis script
- `R/method_palette.R` — shared figure colours
- `CB2/` — pinned submodule, used only as a benchmark baseline and data source
- `examples/` — quickstart, benchmark scripts, and manuscript figures
- `scripts/`, `data/derived/`, `results/` — data preparation and versioned metrics
- `overleaf/` — manuscript LaTeX submodule, mirroring the Overleaf project
- `figures/` — generated figures; the three the manuscript cites are published
  into `overleaf/figures/` by `scripts/publish_figures.sh`

## Reproduce

```sh
Rscript -e 'devtools::test("BARCS")'
Rscript examples/barcs_quickstart.R
latexmk -pdf -cd overleaf/main.tex
```

The package's own test suite lives with the package, in
`BARCS/tests/testthat/`. It replaced the former `tests/run_tests.R` in this
repository when BARCS was split out.

Each benchmark is a single script under `examples/` that writes versioned
outputs to `results/` or `data/derived/`. Some need external tools: official
MAGeCK 0.5.9.5 at `.venv/bin/mageck`, and Julia for the CRISPulator simulation
(`julia/simulate_crispulator_facs.jl`).

## Documentation

- [`docs/barcs-input-output-examples.md`](docs/barcs-input-output-examples.md) —
  input table, function call, output table
- [`docs/barcs-gene-methods.md`](docs/barcs-gene-methods.md) — the four
  guide-to-gene statistics and their diagnostics
- [`docs/barcs-external-method-comparison.md`](docs/barcs-external-method-comparison.md) —
  comparison with MAGeCK-MLE, edgeR-QL, DESeq2, and limma-voom
- [`DEVELOPMENT.md`](DEVELOPMENT.md) — submodule workflow, Overleaf sync, and
  the continuous PDF build

## Scope and limitations

Use CB² for a direct two-condition comparison, and a specialist joint model
such as Waterbear when several bins are correlated partitions of the same
biological pool. BARCS is a prototype for designs where each independently
sequenced library carries a quantitative or multivariable sample-level design.

The covariance is model based and treats the guide-wise dispersion estimate as
a fixed plug-in, so small-sample calibration is not guaranteed: in the
committed continuous-dose simulation, unmoderated BARCS has null type-I error
0.081 and empirical FDP 0.130 at nominal 0.05. Control-tail calibration
requires controls that share one design and must be evaluated on held-out or
cross-fitted controls rather than the guides used to estimate the scale.

BARCS and CB² share a library-total-conditional beta-binomial variance
principle, but their estimators, effect scales, dispersion assumptions, and
reference degrees of freedom differ; no formal equivalence is claimed. The
benchmarks support further evaluation of the moderated and control-denominator
variants. They do not establish superiority over negative-binomial or general
RNA-sequencing methods, and the guide hierarchy does not turn multiple
reagents into biological replicates.
