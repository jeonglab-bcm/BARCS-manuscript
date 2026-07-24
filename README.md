# CB²-Reg: beta-binomial regression for quantitative CRISPR screens

CB² answers a two-condition question: does guide abundance differ between two
groups? CB²-Reg answers the more general coefficient question: does guide
abundance follow dose, time, an ordered phenotype, or an adjusted treatment
effect?

The extension preserves CB²'s central sampling idea—binomial sequencing
variation plus between-library beta-binomial heterogeneity—but replaces the
two-group mean with an ordinary R model matrix. Guide counts remain the
response. Continuous phenotypes, donors, batches, interactions, and other
sample variables are predictors. Coefficients and named contrasts use a
Student t reference based on independent samples rather than reads.

This is an additive path inside CB2, not a replacement for its existing
two-group workflow:

| Feature | Original CB² | CB²-Reg |
|---|---|---|
| Primary question | Difference between two groups | Trend, adjusted association, interaction, or contrast |
| Design | Reference/comparison labels | Arbitrary full-rank model matrix |
| Effect | Group difference or fold change | Conditional log-odds coefficient |
| Variance | Beta-binomial | Same beta-binomial principle |
| Adjustment | Limited by pairwise design | Batch, donor, dose, time, splines, interactions |
| API | Existing CB2 functions | `bbreg()`, `bb_contrast()`, `bb_screen()` |

Official MAGeCK is kept external as a real negative-binomial benchmark; this
repository does not embed a partial reimplementation.

## Repository architecture

This is the umbrella analysis repository. The R package lives in
[`jeonglab-bcm/CB2`](https://github.com/jeonglab-bcm/CB2) and is tracked here
as the `CB2/` Git submodule. The parent repository records the exact CB2 commit
used for every manuscript and benchmark revision.

Clone both histories together:

```sh
git clone --recurse-submodules https://github.com/jeonglab-bcm/CB2-Reg.git
cd CB2-Reg
```

Package changes are committed and pushed from inside `CB2/`; the parent then
commits the updated submodule pointer. See
[`DEVELOPMENT.md`](DEVELOPMENT.md) for the complete two-repository workflow.

## Contents

- `R/bbreg.R`: dependency-free R implementation for one guide, contrasts, and
  guide-by-guide screens.
- `examples/simulation.R`: reproducible beta-binomial simulation comparing the
  beta-binomial t test, a misspecified binomial z test, and official MAGeCK-MLE.
- `examples/gse70038_comparison.R`: head-to-head analysis of all 64,747 guides
  in GSE70038 using a Table-5-style design and official MAGeCK-MLE 0.5.9.5.
- `examples/sanson_benchmark.R`: independent essential-versus-nonessential-gene
  benchmark on the Sanson A375 Brunello screen bundled with CB2, with and
  without official MAGeCK piecewise CNV correction.
- `examples/chronos_tzelepis_benchmark.R`: genuinely longitudinal HT-29
  benchmark (pDNA plus days 7, 10, 13, 16, 19, 22, and 25) against an official
  continuous-time MAGeCK-MLE fit and the deposited Chronos, MAGeCK, and BAGEL2
  results.
- `examples/waterbear_facs_benchmark.R`: ordered four-bin GSE242880 IL2RA
  FACS stress test against official all-bin MAGeCK-MLE, paper-matched
  outer-bin MAGeCK, and the published Waterbear and MAUDE validation results.
- `data/derived/A375_DepMap19Q3_CNV.tsv`: gene-level A375 copy-number profile
  extracted from DepMap Public 19Q3 (ACH-000219).
- `data/derived/HT29_DepMap20Q2_CNV.tsv`: gene-level HT-29 copy-number profile
  extracted from DepMap Public 20Q2 (ACH-000552).
- `scripts/mageck_compat.py` and `scripts/mageck_cnv_correct.py`: runtime
  compatibility and direct access to MAGeCK 0.5.9.5's official CNV normalizer.
- `CB2/`: clone of `jeonglab-bcm/CB2` with the additive `bbreg()`,
  `bb_contrast()`, `bb_screen()`, and negative-control calibration API,
  RcppArmadillo weighted-IRLS kernels, package tests, documentation, and a
  continuous-phenotype vignette.
- `tests/run_tests.R`: base-R regression and input-validation tests.
- `main.tex`: derivation, interpretation, simulation results, and limitations.
- `output/pdf/beta-binomial-regression-continuous-phenotypes.pdf`: rendered
  manuscript.

## Reproduce

From the repository root:

```sh
Rscript tests/run_tests.R
Rscript examples/simulation.R
Rscript examples/gse70038_comparison.R
Rscript examples/sanson_benchmark.R
Rscript examples/chronos_tzelepis_benchmark.R
Rscript examples/waterbear_facs_benchmark.R
Rscript -e 'devtools::test("CB2")'
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The GSE70038 script expects official MAGeCK 0.5.9.5 at
`.venv/bin/mageck`. A project-local installation can be prepared with:

```sh
python3 -m venv .venv
.venv/bin/pip install numpy scipy
curl -L -o /tmp/mageck-0.5.9.5.tar.gz \
  'https://sourceforge.net/projects/mageck/files/0.5/mageck-0.5.9.5.tar.gz/download'
.venv/bin/pip install /tmp/mageck-0.5.9.5.tar.gz
```

The simulation writes its numeric results to `results/`, its plot to `figures/`,
and prints the key operating characteristics. With the included seed, the
beta-binomial t test has null type-I error 0.081, power 1.000, and empirical
FDR 0.130. Official MAGeCK-MLE has type-I error 0.038, power 0.975, and
empirical FDR 0.049. The read-level binomial z test has type-I error 0.844 and
empirical FDR 0.770. Table 1 and Figure 1 in the manuscript show all three
methods.

On GSE70038, beta-binomial versus official MAGeCK-MLE gene-effect Spearman
correlations are 0.888-0.928 across the four terminal-condition coefficients;
top-200 depleted-gene Jaccard overlaps are 0.533-0.562. The analysis uses
MAGeCK's Wald p-values/FDR and writes all guide, gene, concordance, and
published-validation-gene tables under `results/gse70038/`.
The rank correlation is computed here, not taken from GEO or a paper: for each
coefficient, within-gene median beta-binomial guide effects are matched to
official MAGeCK gene betas, both vectors are ranked, and the paired ranks are
correlated. Every contributing pair is in
`results/gse70038/effect_concordance_pairs.csv.gz`. It measures agreement, not
which model fits better.

On the independent Sanson A375 gold-standard benchmark, beta-binomial and
official MAGeCK-MLE have essentially identical uncorrected AUROC (0.9598 on
the CNV-complete gene set). MAGeCK's official piecewise CNV adjustment was
applied to both methods' effects. Among reference nonessential genes, the
effect–CNV Spearman correlation changes from -0.185 to -0.042 for
beta-binomial and from -0.207 to -0.006 for MAGeCK, showing successful
removal of CNV-associated depletion. CNV correction does not improve global
essential-gene AUROC in this dataset and does not change fixed-FDR calls
because MAGeCK 0.5.9.5 adjusts beta scores after calculating p-values and FDR.

After CNV correction, beta-binomial retains a recall advantage of 0.0418 and
an F1 advantage of 0.0253 at nominal FDR 0.05; both paired-bootstrap intervals
exclude zero. Reproducible scores, CNV diagnostics, threshold curves, metrics,
and bootstrap intervals are under `results/sanson_benchmark/`.

The Tzelepis/Chronos benchmark is the direct test of the continuous-time use
case. The three sequencing columns at each day are summed before inference, so
technical replication does not inflate the t-test degrees of freedom. On the
shared set of 1,080 reference-essential and 5,774 unexpressed HT-29 genes,
beta-binomial time regression, official continuous-time MAGeCK-MLE, and
Chronos-joint have AUROC 0.9786, 0.9789, and 0.9747, respectively. MAGeCK has
the highest PR AUC (0.9534), beta-binomial has the highest recall subject to
at least 90% precision (0.9037), and Chronos has the strongest normalized null-median
difference (-18.88). Thus no method wins every target: Chronos best separates
the distribution centers, MAGeCK slightly leads global ranking, and the
beta-binomial slope leads the high-precision recall criterion. The
beta-binomial and official MAGeCK numeric-time effects have Spearman
correlation 0.926; as in the GSE70038 analysis, this is computed agreement, not
a likelihood comparison.

This dataset also exposes a limitation of post-fit single-line CNV
normalization. The unexpressed-gene effect–CNV correlation changes from -0.072
to -0.145 for beta-binomial and from -0.075 to -0.149 for MAGeCK after applying
MAGeCK 0.5.9.5's official piecewise correction. In other words, the correction
over-adjusts this longitudinal screen instead of removing the already weak
association. Chronos was left uncorrected because its published CNV procedure
requires multiple cell lines. The source effects, common evaluation universe,
metrics, rank correlations, and CNV audit are under
`results/chronos_tzelepis/`.

The GSE242880 benchmark uses the low-coverage, high-MOI primary-T-cell arm:
four ordered IL2RA FACS bins for each of three donors. The 26 evaluation genes
were validated by individual knockout and flow cytometry. Raw all-bin
CB²-Reg recovers 23/26 at gene FDR 0.10 but calls 127 genes
and is inflated among 593 non-targeting guides (13.3% have nominal
guide-level p < 0.05). The new `bb_calibrate_controls()` tail-scale diagnostic
restores that frequency to 4.9%; the calibrated result recovers 22/26 with 49
discoveries. Official all-bin MAGeCK-MLE recovers 17/26 with 72 discoveries,
and the independent paper-matched outer-bin `mageck test` rerun recovers 18/26
with 32. Waterbear's published, not rerun, result is 24/26 with 79 calls;
MAUDE's is 25/26 with 406 calls.

The deposited follow-up table also contains seven candidates that did not
validate experimentally. In this selected 33-gene panel, calibrated all-bin
CB²-Reg has F1 0.863, Matthews correlation 0.398, balanced accuracy 0.709,
AUROC 0.808, and average precision 0.945. All-bin MAGeCK-MLE has 0.723, 0.070,
0.541, 0.626, and 0.872, respectively. These are supporting metrics from a
small candidate-selected panel, not unbiased genome-wide accuracy estimates.
They nevertheless show why the calibrated result is preferable to raw CB²-Reg:
the continuous-bin signal is retained while the null tail is repaired.

FACS bins are correlated partitions. Waterbear models that joint structure and
remains the best-fitting method for this experimental design; CB²-Reg provides
a fast, transparent trend analysis and a useful sensitivity check. Full tables
are under `results/waterbear_facs/`, including
`validation_panel_metrics.csv`.

The CB2 package benchmark processes roughly 2,100 guides/second serially and
7,375 guides/second with four forked workers on the current machine (about
1.36 seconds for 10,000 guides, or 8.8 seconds projected for 64,747 guides).
The new full longitudinal fit processes 86,882 filtered guides in about 20
seconds with four workers on the same machine.

## Important scope note

Use original CB² for a direct two-condition comparison. Use CB²-Reg when each
independently sequenced library has a quantitative or multivariable
sample-level design. Use a specialist joint model such as Waterbear when
several bins are correlated partitions of the same biological pool.

CB²-Reg is not a model for a continuous phenotype measured per cell when guide
identity and phenotype are not jointly observed. The current regression layer
is a transparent research implementation; the manuscript lists the
dispersion-shrinkage, gene-level hierarchy, and repeated-measure work still
needed for broader production use.
