# BARCS Web

BARCS Web is a static, privacy-first implementation of BARCS. The count
matrix, model fitting, calibration, plots, and result export all remain in the
browser. It does not require an R or Python installation, a database, or an
analysis server.

## What it implements

- comma- or tab-separated guide-count and sample-metadata input;
- full-library totals supplied through a metadata `total`, `library_total`, or
  `library_size` column, or computed from all uploaded guide rows;
- numeric and treatment-coded categorical predictors;
- additive covariates and predictor-by-covariate interactions;
- guide-wise beta-binomial feasible IRLS;
- Student-t coefficient tests using sample residual degrees of freedom;
- Benjamini-Hochberg FDR;
- optional negative-control tail calibration;
- optional shared-effect gene statistics without Fisher or Stouffer
  aggregation;
- convergence, dispersion, denominator, and calibration diagnostics;
- downloadable guide- and gene-level CSV files.

The numerical engine is in `public/barcs-core.js` and runs in
`public/barcs-worker.js`, keeping large analyses off the interface thread.

## Input contract

The count table must contain `guide` and one column for every sample. `gene`
and `control` are optional but required for gene inference and control
calibration, respectively.

```csv
guide,gene,control,D0_A,D0_B,D7_A,D7_B
g1,GENE1,false,120,98,54,47
nt1,NTC,true,85,91,87,90
```

The metadata table must contain `sample`. Other columns are available as
predictors. Supplying the optional full-library `total` is recommended when
the uploaded guide table has already been filtered.

```csv
sample,time,batch,total
D0_A,0,A,1824000
D0_B,0,B,1912500
D7_A,7,A,1768300
D7_B,7,B,1897100
```

Sample names must match the count columns exactly. The number of independent
libraries must exceed the number of fitted design columns.

## Develop and verify

```sh
cd web
npm run dev
npm test
```

The test suite checks parsing, design coding, probability calculations,
calibration, gene inference, the deployment artifact, and numerical parity
with reference values produced by `R/bbreg.R`. The three parity fixtures agree
with the R coefficient, standard error, t statistic, p-value, and dispersion
within numerical tolerance.

## Deployment

`npm run build` produces a dependency-free static client under `dist/client`
and a Cloudflare-compatible worker entry under `dist/server/index.js`. No
durable storage or server-side analysis bindings are declared.

## Scope

The web implementation reproduces the current BARCS guide coefficient model.
It does not turn sequencing reads or multiple guides into biological
replicates, reconstruct unobserved single-cell phenotypes, or model dependence
among FACS bins derived from the same biological pool.
