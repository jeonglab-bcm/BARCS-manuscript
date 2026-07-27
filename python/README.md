# BARCS in Python

`barcs` is an rpy2 wrapper around [`R/bbreg.R`](../R/bbreg.R). It reimplements
nothing: every function calls into that file.

```sh
pixi run test-python   # marshalling checks + a full run against native Rscript
pixi run notebook      # interactive walkthrough
```

```python
import sys; sys.path.insert(0, "python")
from barcs import bbreg, bb_screen, bb_gene_eb_moderate

fit = bbreg(counts_for_one_guide, library_totals, "~ dose + batch", samples)
print(fit.summary())

screen = bb_screen(counts, samples, "~ dose + batch", term="dose",
                   totals=library_totals, guide=guides, gene=genes)
genes = bb_gene_eb_moderate(screen)
```

## Why a wrapper and not a port

An earlier version of this module was a from-scratch NumPy/SciPy
reimplementation, held to the original by a test suite requiring the two to
agree to eight significant figures. It worked, and it was still the wrong
shape: a second place for the statistics to live and drift, covering only five
of the ten public functions — the ones judged worth reimplementing.

Calling R directly means there is exactly one implementation, and **all ten
functions are available**, including the three later guide-to-gene statistics
and `bb_moderate_dispersion` that the manuscript's later comparisons actually
turn on. A change to `R/bbreg.R` is picked up on the next call and cannot
silently disagree with what the benchmarks used.

The cost is an R installation with Rcpp/RcppArmadillo — already repository
dependencies declared in `pixi.toml`. This does not add a requirement so much
as stop pretending the Python side didn't have one.

## What is available

| | |
|---|---|
| `bbreg` | single-guide fit, returned as `BBRegResult` |
| `bb_contrast` | linear contrasts |
| `bb_screen` | guide-by-guide screen with BH-adjusted FDR |
| `bb_calibrate_controls` | empirical-null rescaling, both estimators |
| `bb_moderate_dispersion` | guide-dispersion moderation |
| `bb_gene_original` | signed-z aggregation |
| `bb_gene_normal` | exchangeable normal guide coefficients |
| `bb_gene_partial_pool` | random-effects partial pooling |
| `bb_gene_eb_moderate` | empirical-Bayes heterogeneity moderation |
| `bb_gene_consistency` | guide consistency for single-replicate screens |

## How it is tested

There is nothing statistical to verify — the numerics are R's. What can break
is the marshalling, so that is what the tests pin.

`test_matches_rscript.py` fits one fixture screen twice: once through the
wrapper, once by running `reference_fit.R` under a plain `Rscript` over CSV
files. Both source the same `R/bbreg.R` with no CB2 loaded, so they execute
identical code and **must agree to 1e-12** — CSV round-trip noise and nothing
else. When this module was a port the bar was 1e-7, with the residual coming
from R's `uniroot` and SciPy's `brentq` stopping in different places; there is
now only one root finder.

`test_r_bridge.py` pins the individual conversions that fail quietly: a
guide-by-sample matrix arriving transposed (R fills matrices column-major), a
pandas Categorical losing its level order and silently changing the treatment
coding, R attributes such as `control_scale` being dropped, and each of the
four gene statistics actually being reachable.

## Things worth knowing

**Results are pandas, R objects are still there.** Data frames come back as
DataFrames with a clean `RangeIndex` — R's `"1"`, `"2"`, … row names are
dropped, so `screen["guide"][0]` works rather than raising. R attributes land
in `.attrs`: `bb_calibrate_controls` puts the null scale in
`.attrs["control_scale"]`. `BBRegResult` keeps the untouched R fit on
`.r_object`, which is how `bb_contrast` passes a fit back into R without
refitting.

**`min_total_count` defaults to 10**, R's default, not the old port's 1.

**Formulas only.** `bbreg` and `bb_screen` no longer accept a pre-built design
matrix; pass `"~ dose + batch"` and let R's formula machinery handle it. That
is what keeps coefficient names — and therefore `term=` — identical to a
native R session.

**Importing costs an R session.** `python/barcs/bbreg.py` sources
`R/bbreg.R` at import time, so `import barcs` starts an embedded R.

## The notebook

`notebooks/barcs_walkthrough.py` is a [marimo](https://marimo.io) notebook,
stored as plain Python so it diffs and reviews like source. It is reactive:
moving a slider re-runs exactly the cells that depend on it.

It simulates a dose-response screen with known truth, then walks through a
single fit, a contrast, the full screen, control calibration, and all four
guide-to-gene statistics — alongside a plain binomial GLM on the same design.
At the default settings the binomial rejects about **85%** of true nulls at
p<0.05 where BARCS rejects about 9%; setting `rho` to zero collapses the gap,
which is the check that the difference comes from overdispersion and not from
something else in the comparison.

The notebook is explicit that 9% is not 5%, and why: the covariance treats
each guide's dispersion as a fixed plug-in, so the t reference does not carry
the uncertainty in having estimated it.

It is a demonstration, not a benchmark. The benchmarks are the R scripts in
`examples/`; [`docs/repository-map.md`](../docs/repository-map.md) says which
one answers what.
