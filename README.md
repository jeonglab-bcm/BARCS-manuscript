# BARCS: Beta-binomial Analysis and Regression for CRISPR Screens

CB² answers a two-condition question: does guide abundance differ between two
groups? BARCS (Beta-binomial Analysis and Regression for CRISPR Screens)
answers the more general coefficient question: does guide abundance follow
dose, time, an ordered phenotype, or an adjusted treatment effect?

The extension preserves CB²'s central sampling idea—binomial sequencing
variation plus between-library beta-binomial heterogeneity—but replaces the
two-group mean with an ordinary R model matrix. Guide counts remain the
response. Continuous phenotypes, donors, batches, interactions, and other
sample variables are predictors. Coefficients and named contrasts use a
Student t reference based on independent samples rather than reads.

The mathematical relationship is deliberately narrow. A legacy representation
lemma shows that a saturated two-cell generalized least-squares calculation
reproduces original CB² after its weighted group proportions and variances have
already been computed. This is an algebraic compatibility identity, not proof
that default `bbreg()` strictly nests the original estimator. Under a correctly
specified common-dispersion binary design, local alternatives, and increasing
numbers of independent libraries, the logit-scale BARCS statistic is
first-order equivalent to the raw-proportion statistic. Increasing read depth
at fixed biological replication is not sufficient. See
`examples/cb2_generalization_proof.R` for the machine-precision representation
check and convergence illustration.

The current covariance is model based and treats the guide-wise
dispersion estimate as a fixed plug-in value. Residual Student t degrees of
freedom do not formally propagate dispersion-estimation uncertainty, so
small-sample calibration remains a stated limitation. Control-tail calibration
also requires controls to share one residual degree of freedom, as they do
under one common complete design.

For a version written without matrix algebra, start with
[`docs/cb2-generalization-high-school-proof.md`](docs/cb2-generalization-high-school-proof.md).
For a code-first version with six progressively richer examples, use
[`docs/cb2-generalization-computational-walkthrough.md`](docs/cb2-generalization-computational-walkthrough.md).
For compact input-table → function-call → output-table examples, use
[`docs/barcs-input-output-examples.md`](docs/barcs-input-output-examples.md).
For installation-free analysis, [`web/`](web/) contains BARCS Web: a static
browser implementation that fits the same beta-binomial estimating equations
locally, validates the design, draws diagnostics, and exports guide- and
gene-level results without sending screen data to a server.

This is an additive path inside CB2, not a replacement for its existing
two-group workflow:

| Feature | Original CB² | BARCS |
|---|---|---|
| Primary question | Difference between two groups | Trend, adjusted association, interaction, or contrast |
| Design | Reference/comparison labels | Arbitrary full-rank model matrix |
| Effect | Group difference or fold change | Conditional log-odds coefficient |
| Dispersion | Estimated separately within groups | One guide-wise rho across the fitted design |
| Estimator | Weighted group proportions | Feasible logistic IRLS |
| Degrees of freedom | Welch–Satterthwaite | Residual sample df |
| Variance principle | Beta-binomial | Same beta-binomial principle |
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
git clone --recurse-submodules https://github.com/jeonglab-bcm/BARCS.git
cd BARCS
```

Package changes are committed and pushed from inside `CB2/`; the parent then
commits the updated submodule pointer. See
[`DEVELOPMENT.md`](DEVELOPMENT.md) for the complete two-repository workflow.

## Contents

- `R/bbreg.R`: dependency-free R implementation for one guide, contrasts, and
  guide-by-guide screens.
- `web/`: dependency-free serverless web application with a browser worker,
  local FASTQ/gzip alignment, candidate-library determination, mapping and
  representation QC, downloadable results, and full-example plus
  manuscript-scale parity tests against the reference R implementation. It
  includes a compact Liang HAP1 processed-count example and a matching
  synthetic FASTQ teaching bundle built from real Liang guide sequences. The
  two implementations are scientifically equivalent within declared
  floating-point tolerances, not promised to be bit-for-bit identical.
- `julia/simulate_crispulator_facs.jl`: pinned CRISPulator 0.5.1 simulation of
  low 25%, high 25%, overlapping 0--100% bulk, and input samples.
- `examples/crispulator_facs_benchmark.R`: one-seed BARCS/MAGeCK-MLE FACS
  analysis with truth-based ranking, directional recovery, and null metrics.
- `examples/crispulator_facs_repeated_benchmark.R`: five-seed manuscript
  benchmark and aggregate figure. The older `examples/simulation.R` remains as
  a beta-binomial calibration diagnostic but is no longer a main result.
- `examples/cb2_generalization_proof.R`: numerical verification of the legacy
  GLS representation and a local-equivalence illustration for the logit
  statistic.
- `examples/cb2_generalization_walkthrough.R`: six code-first examples,
  including hand-checkable numbers, 1,000 randomized identity checks, original
  `fit_ab()` counts, the finite-sample logit difference, local convergence, and
  a continuous-dose fit.
- `examples/barcs_input_output_examples.R`: printable two-group,
  continuous-dose, and dose-plus-batch inputs with their exact model outputs.
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
- `examples/liang_cas13_benchmark.R`: five-cell-line Cas13 fitness
  processed-count sensitivity analysis using Liang's deposited
  RobustRankAggreg results, official MAGeCK-RRA, official MAGeCK-MLE, and
  BARCS on the same normalized day-0/day-14 values.
- `scripts/prepare_liang_cas13.R`,
  `scripts/count_liang_cas13_run.sh`, and
  `scripts/queue_liang_cas13_counts.sh`: download the Liang supplementary
  tables, stream the 20 endpoint FASTQs through the published anchor/Bowtie
  rules without retaining reads, and submit restartable counting jobs.
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
Rscript examples/cb2_generalization_proof.R
Rscript examples/cb2_generalization_walkthrough.R
Rscript examples/barcs_input_output_examples.R
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
julia --project=julia julia/simulate_crispulator_facs.jl
Rscript examples/crispulator_facs_repeated_benchmark.R
Rscript examples/gse70038_comparison.R
Rscript examples/sanson_benchmark.R
Rscript examples/chronos_tzelepis_benchmark.R
Rscript examples/waterbear_facs_benchmark.R
Rscript scripts/prepare_liang_cas13.R
Rscript examples/liang_cas13_benchmark.R
Rscript -e 'devtools::test("CB2")'
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The primary Liang comparison deliberately uses the deposited normalized
guide-count values, as they make the processed-data comparison quick and fully
reproducible. Liang's published “RRA” is `RobustRankAggreg` 1.2.1, not
MAGeCK-RRA, so the benchmark keeps those as separate methods. The deposited
values are fractional after median-ratio normalization, ComBat correction,
and outlier processing. Because BARCS enforces integer counts, the script
explicitly rounds them to nearest pseudo-counts and gives that identical
matrix to the three newly fitted methods--BARCS, MAGeCK-RRA, and MAGeCK-MLE;
it also records the rounding error. Consequently, this is labelled a
processed-count sensitivity analysis:
neither BARCS nor MAGeCK has its literal raw-count sampling likelihood.
Known-essential protein-coding controls are positives
and cell-line-specific non-expressed lncRNAs are null controls. Published RRA
calls are a comparator, not a circular definition of truth.

The processed-count result is intentionally not a BARCS win. Macro-averaged
over the five cell lines, Liang RRA has the highest AUROC (0.970), average
precision (0.891), essential recall at 5% null FPR (0.917), and directional
essential recall at FDR 0.10 (0.823). BARCS reaches 0.939, 0.776, 0.823, and
0.487, respectively. Every BARCS day-effect fit has only one residual degree
of freedom after the available replicate block, so its main limitation is
threshold power rather than a complete loss of biological ranking. Mean
BARCS--MAGeCK-MLE effect-rank correlation is 0.875 across cell lines. The
versioned metrics and rounding audit are under `data/derived/`.

An optional raw-read confirmation remains available. It streams approximately
14 GB of compressed endpoints without retaining FASTQs:

```sh
bash scripts/queue_liang_cas13_counts.sh 2
pueue wait --group liang-cas13
```

CRISPulator MOI, guide quality, library size, and replicate count are
configurable. Quality is the fraction of high-quality CRISPRn guides:
high-quality guides have complete knockout in the simulator, whereas
low-quality guide activity follows CRISPulator's truncated low-activity
distribution. The manuscript baseline is MOI 0.25, 90% high-quality guides,
400 genes, and four independent screen replicates.

```sh
# The default sensitivity analysis is ten scenarios x five seeds:
# three MOIs, guide-quality fractions, gene counts, and replicate counts are
# varied one at a time (50 simulations, including the baseline and a
# diagnostic one-replicate boundary case).
Rscript examples/crispulator_facs_repeated_benchmark.R

# Optional 3 x 3 x 3 x 4 factorial (540 simulations across five seeds).
CRISPULATOR_GRID_MODE=full_factorial \
CRISPULATOR_MOI_VALUES=0.10,0.25,0.40 \
CRISPULATOR_HIGH_QUALITY_GUIDE_FRACTION_VALUES=0.60,0.75,0.90 \
CRISPULATOR_GENE_VALUES=100,400,1000 \
CRISPULATOR_REPLICATE_VALUES=1,3,4,6 \
Rscript examples/crispulator_facs_repeated_benchmark.R

# Run one custom scenario across five seeds.
CRISPULATOR_GRID_MODE=single \
CRISPULATOR_MOI=0.40 \
CRISPULATOR_HIGH_QUALITY_GUIDE_FRACTION=0.60 \
CRISPULATOR_GENES=1000 \
CRISPULATOR_REPLICATES=6 \
Rscript examples/crispulator_facs_repeated_benchmark.R

# Equivalent direct Julia arguments:
# output directory, replicates, seed, MOI, high-quality-guide fraction, genes
julia --project=julia julia/simulate_crispulator_facs.jl \
  results/custom_facs 4 20250729 0.40 0.60 1000
```

The multimethod benchmark requires the Bioconductor packages `edgeR`,
`DESeq2`, and `limma` in addition to the official MAGeCK executable:

```r
BiocManager::install(c("edgeR", "DESeq2", "limma"))
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

The CRISPulator baseline uses five fixed seeds, 400 genes, five guides per
gene, and four independent screen replicates. In the low--bulk--high design,
MAGeCK-MLE slightly leads global average precision (0.924 versus 0.917), while
BARCS leads directionally correct recall at gene FDR 0.10 (0.764 versus 0.717)
and F1 (0.828 versus 0.817). Low--bulk--high and tail-only effect ranks are
nearly identical (mean Spearman 0.9985 for BARCS and 0.9999 for MAGeCK-MLE).
Bulk versus input remains null (mean AUROC 0.500 and effect Spearman 0.060).
Thus bulk can stabilize variance and degrees of freedom but adds no directional
FACS contrast; it also overlaps the tails and uses 50% more sequencing.
The same continuous design is fitted with BARCS, MAGeCK-MLE, edgeR-QL,
DESeq2, and limma-voom. The latter three use a common directional Stouffer
guide-to-gene summary, whereas MAGeCK-MLE retains its native gene model.
Across the nine supported scenarios with at least three replicates,
BARCS-minus-MAGeCK average
precision ranges from -0.015 to 0.022, showing no universal ranking advantage.
The directional-recall difference is positive in seven of nine scenarios and
ranges from -0.018 to 0.136; the F1 difference ranges from -0.022 to 0.080.
Across those nine scenarios, realized FDP averages 0.094 for BARCS and 0.064
for MAGeCK-MLE, versus 0.230-0.250 for the three general count-model pipelines.
The additional one-replicate stress test is reported separately. Ordinary
BARCS improves average precision over MAGeCK-MLE (0.633 versus 0.548) but
makes no discoveries at gene FDR 0.10. The optional `bb_gene_consistency()`
analysis (`BARCS-GC`) estimates a shared guide coefficient by inverse-variance
weighting and calibrates its Wald statistic with a robust gene-level empirical
null; it does not use Fisher or Stouffer aggregation. Across the same five
seeds it raises average precision to 0.670, directional recall to 0.145, and
F1 to 0.245, with mean realized FDP 0 in the simulation. This remains a
hypothesis-ranking tool: multiple guides demonstrate perturbation
reproducibility but do not replace biological replication.

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
BARCS recovers 23/26 at gene FDR 0.10 but calls 127 genes
and is inflated among 593 non-targeting guides (13.3% have nominal
guide-level p < 0.05). The new `bb_calibrate_controls()` tail-scale diagnostic
restores that frequency to 4.9%; the calibrated result recovers 22/26 with 49
discoveries. Official all-bin MAGeCK-MLE recovers 17/26 with 72 discoveries,
and the independent paper-matched outer-bin `mageck test` rerun recovers 18/26
with 32. Waterbear's published, not rerun, result is 24/26 with 79 calls;
MAUDE's is 25/26 with 406 calls.

The deposited follow-up table also contains seven candidates that did not
validate experimentally. In this selected 33-gene panel, calibrated all-bin
BARCS has F1 0.863, Matthews correlation 0.398, balanced accuracy 0.709,
AUROC 0.808, and average precision 0.945. All-bin MAGeCK-MLE has 0.723, 0.070,
0.541, 0.626, and 0.872, respectively. These are supporting metrics from a
small candidate-selected panel, not unbiased genome-wide accuracy estimates.
They nevertheless show why the calibrated result is preferable to raw BARCS:
the continuous-bin signal is retained while the null tail is repaired.

FACS bins are correlated partitions. Waterbear models that joint structure and
remains the best-fitting method for this experimental design; BARCS provides
a fast, transparent trend analysis and a useful sensitivity check. Full tables
are under `results/waterbear_facs/`, including
`validation_panel_metrics.csv`.

The CB2 package benchmark processes roughly 2,100 guides/second serially and
7,375 guides/second with four forked workers on the current machine (about
1.36 seconds for 10,000 guides, or 8.8 seconds projected for 64,747 guides).
The new full longitudinal fit processes 86,882 filtered guides in about 20
seconds with four workers on the same machine.

## Important scope note

Use original CB² for a direct two-condition comparison. Use BARCS when each
independently sequenced library has a quantitative or multivariable
sample-level design. Use a specialist joint model such as Waterbear when
several bins are correlated partitions of the same biological pool.

BARCS is not a model for a continuous phenotype measured per cell when guide
identity and phenotype are not jointly observed. The current regression layer
is a transparent research implementation; the manuscript lists the
dispersion-shrinkage, gene-level hierarchy, and repeated-measure work still
needed for broader production use.
