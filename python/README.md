# BARCS in Python

A Python implementation of the BARCS guide-level regression, plus a marimo
notebook that drives it interactively.

```sh
pixi run test-python-vs-r   # check it against R
pixi run notebook           # interactive walkthrough
```

```python
import sys; sys.path.insert(0, "python")
from barcs import bbreg, bb_screen, bb_gene_original

fit = bbreg(counts_for_one_guide, library_totals, "~ dose + batch", samples)
print(fit.summary())

screen = bb_screen(counts, samples, "~ dose + batch", term="dose",
                   totals=library_totals, guide=guides, gene=genes)
genes = bb_gene_original(screen)
```

## What is here, and what is not

| | |
|---|---|
| `bbreg` | single-guide fit, full coefficient table |
| `bb_contrast` | linear contrasts, with the covariance carried through |
| `bb_screen` | guide-by-guide screen with BH-adjusted FDR |
| `bb_calibrate_controls` | empirical-null rescaling, both estimators |
| `bb_gene_original` | the signed-z guide-to-gene statistic |

**Not ported**: the exchangeable-normal, partial-pooling, and empirical-Bayes
guide-to-gene statistics, and `bb_moderate_dispersion`. Those are what the
manuscript's later comparisons turn on, and every committed result involving
them came from R. This is a deliberate boundary, not an unfinished edge — use
`R/bbreg.R` for that work rather than assuming the Python module covers it.

## Why a port rather than a wrapper

Calling R from Python through rpy2 would have kept a single implementation, at
the cost of making every Python user install R, the CB2 submodule, and a
compiler. This module has no R dependency at all.

That only means something if the two agree, so `tests/test_against_r.py`
builds one fixture screen, fits it in both languages, and requires agreement
across every reported quantity. It is not a smoke test: it compares each
coefficient of a single fit, a contrast, all 96 guides of a screen including
the BH-adjusted FDR, both control-calibration estimators, and the gene-level
aggregation.

## How closely they agree, and why not exactly

Coefficient estimates match to about **1e-12**, which is double-precision
noise. Everything downstream is slightly looser, for one traceable reason.

Both implementations estimate `rho` by solving `pearson(rho) = df_residual`.
R uses `uniroot(tol = 1e-10)`, which on the fixture stops where the Pearson
statistic is `11.0000000083` against a target of `11`. SciPy's `brentq`, given
a tighter tolerance, reaches `10.999999999999776`. The two therefore land on
`rho` values differing by 7.5e-10 — and **the Python root is the more accurate
one**.

That difference propagates, growing a little wherever one noisy quantity is
divided by another:

| quantity | agreement |
|---|---|
| coefficient estimates | 1e-12 |
| standard errors | 4e-10 |
| negative-control scale | 1e-8 |
| calibrated FDR | 3e-8 |

The test threshold is 1e-7 — eight significant figures on every reported
quantity. It is not tighter because tightening it would mean deliberately
degrading the Python solver to imitate R's stopping rule.

## Two behaviours worth knowing before you rely on this

**Coefficient names follow R, not formulaic.** `formulaic` writes `Intercept`
and `batch[T.b2]`; this module renames them to R's `(Intercept)` and
`batchb2`. Without that, the same analysis would need a different `term=`
argument in each language and every cross-check against an R result would have
to translate. Pass a design matrix instead of a formula string if you want to
skip formula parsing entirely.

**`converged` is not an abundance filter.** On a guide with only a handful of
reads, R declines to fit at all while Python often converges and returns a
finite estimate. Neither number means anything at that depth. Filter on
abundance — `min_total_count`, or `mean_cpm` afterwards — and do not compare
the two implementations on guides like these.

## The notebook

`notebooks/barcs_walkthrough.py` is a [marimo](https://marimo.io) notebook,
stored as plain Python so it diffs and reviews like source rather than like a
`.ipynb`. It is reactive: moving a slider re-runs exactly the cells that depend
on it, so there is no stale-state problem.

It simulates a dose-response screen with known truth and walks through a single
fit, a contrast, the full screen, and control calibration — alongside a plain
binomial GLM on the same design. At the default settings the binomial rejects
about **85%** of true nulls at p<0.05 where BARCS rejects about 9%; setting
`rho` to zero collapses the gap, which is the check that the difference comes
from overdispersion rather than from something else in the comparison.

It is a demonstration, not a benchmark. The benchmarks are the R scripts in
`examples/`; [`docs/repository-map.md`](../docs/repository-map.md) says which
one answers what.
