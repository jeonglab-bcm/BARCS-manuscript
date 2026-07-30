# BARCS: Beta-binomial Analysis and Regression for CRISPR Screens

CB² asks whether guide abundance differs between two groups. BARCS answers the
general coefficient question: does guide abundance follow dose, time, an
ordered phenotype, or an adjusted treatment effect?

BARCS keeps CB²'s sampling idea—binomial sequencing variation plus
between-library beta-binomial heterogeneity—but replaces the two-group mean
with an ordinary R model matrix. Guide counts remain the response; dose, time,
donor, batch, and interaction terms are predictors. Coefficients and named
contrasts use a Student t reference based on independent samples rather than
reads. This is an additive path inside CB2, not a replacement for its
two-group workflow.

## Install

The R package lives in [`jeonglab-bcm/CB2`](https://github.com/jeonglab-bcm/CB2)
and is tracked here as the `CB2/` submodule, so this repository pins the exact
package commit used for each manuscript and benchmark revision.

```sh
git clone --recurse-submodules https://github.com/jeonglab-bcm/BARCS.git
cd BARCS
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

- `R/bbreg.R` — dependency-free implementation; `R/README.md` is the source map
- `CB2/` — pinned package submodule with the compiled weighted-IRLS kernels
- `examples/` — quickstart, benchmark scripts, and manuscript figures
- `scripts/`, `data/derived/`, `results/` — data preparation and versioned metrics
- `main.tex`, `sections/` — manuscript source; rendered PDF in `output/pdf/`

## Reproduce

```sh
Rscript tests/run_tests.R
Rscript examples/barcs_quickstart.R
Rscript -e 'devtools::test("CB2")'
```

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
- [`DEVELOPMENT.md`](DEVELOPMENT.md) — two-repository submodule workflow

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
