# BARCS Web

BARCS Web is a static, privacy-first FASTQ-to-results implementation of BARCS.
Raw reads, guide libraries, count matrices, model fitting, plots, and result
export all remain in the browser. FASTQ quantification and guide fitting run
in separate web workers. It does not require R, Python, a database, or an
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

## Raw-read workflow

The FASTQ path accepts one or more candidate guide libraries plus one FASTQ or
FASTQ.gz file per independent sample. It:

- parses FASTA or delimited `guide`, `sequence`, `gene`, and `control`
  annotations;
- excludes repeated guide sequences that cannot be assigned uniquely;
- samples the first FASTQ to rank candidate libraries by exact unique matches;
- searches every read in both orientations and determines the dominant guide
  position;
- streams gzip decompression and counting without uploading raw reads;
- reports total, uniquely mapped, ambiguous, and unmapped reads, zero-guide
  fraction, count Gini coefficient, orientation, and guide position;
- produces downloadable counts, QC, and a sample-metadata template.

Exact uniquely mapped guide reads define the full-library totals supplied to
BARCS. FASTQ basenames become sample identifiers and must match the `sample`
column in metadata.

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
with reference values produced by `R/bbreg.R`. The committed suite covers all
48 guides under four supported design shapes (192 fits) and all nine
shared-effect gene results in the example screen. It requires the same guide
and gene order, convergence status, and FDR < 0.10 calls, while comparing
reported statistics under explicit numerical tolerances.

The R implementation is the reference. The browser and R results are
numerically equivalent, not bit-for-bit identical: independent IRLS stopping
paths, matrix arithmetic, root finding, and probability-library
implementations can produce small IEEE-754 floating-point differences. See
[`PARITY.md`](PARITY.md) for the tested contract, observed errors, and
reproduction commands.

## Manuscript verification

The **Verify manuscript** action loads the actual 64,747-guide GSE70038 count
matrix and the 16-by-5 Table 5 design used in the manuscript. All four terminal
coefficients are supported in the same joint model. The preset initially shows
`GSC0131_end`; changing the tested coefficient applies the corresponding
locked reference. It uses the manuscript's explicitly labeled
median-effect/directional-Stouffer comparison summary. The default BARCS gene
result remains the shared-effect Wald construction; Stouffer aggregation is
not required for guide-level BARCS.

For every coefficient, verification checks 64,747 guide rows, 18,077 gene
rows, the complete FDR < 0.10 discovery count, aggregate
effect/p-value/FDR checksums, the top-20 ordering, and six prespecified
validation genes. Across all four coefficients, the largest observed
R-versus-browser differences are `1.40e-11` for guide effects, `8.42e-12` for
guide p-values/FDR, `5.14e-9` for gene p-values, and `5.26e-9` for gene FDR.
No FDR < 0.10 decision differs.

## Liang HAP1 examples

**Load Liang HAP1** opens a compact 72-guide view of the deposited HAP1
processed counts: four guides each for four essential/protein-coding genes,
two published lncRNA signals, and two TPM-zero lncRNA nulls, plus 40
non-targeting controls. The values are rounded normalized/ComBat
pseudo-counts, while metadata totals are retained from the complete
56,322-guide table before subsetting.

**Run Liang FASTQ demo** aligns four bundled FASTQ.gz fixtures against the real
selected Liang guide sequences and then configures `~ replicate + day14`.
These reads are synthetic exact-match teaching fixtures, not deposited Liang
raw reads. They reproduce the bundled selected-guide pseudo-counts exactly and
exercise forward/reverse orientation detection. The complete downloadable
bundle and provenance table are under `public/examples/liang-hap1/`; regenerate
them with:

```sh
Rscript scripts/generate_liang_web_examples.R
```

## Deployment

`npm run build` produces a dependency-free static client under `dist/client`
and a Cloudflare-compatible worker entry under `dist/server/index.js`. No
durable storage or server-side analysis bindings are declared.

## Scope

The exact matcher is intended for pooled-screen amplicons containing the
library spacer. It is not a genomic aligner, adapter trimmer, base-quality
recalibrator, or approximate mismatch caller. The web implementation does not
turn sequencing reads or multiple guides into biological replicates,
reconstruct unobserved single-cell phenotypes, or model dependence among FACS
bins derived from the same biological pool.
