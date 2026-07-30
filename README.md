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
also requires controls to share the same residual degrees of freedom, as they
do under one common complete design.

For a version written without matrix algebra, start with
[`docs/cb2-generalization-high-school-proof.md`](docs/cb2-generalization-high-school-proof.md).
For a code-first version with six progressively richer examples, use
[`docs/cb2-generalization-computational-walkthrough.md`](docs/cb2-generalization-computational-walkthrough.md).
For compact input-table → function-call → output-table examples, use
[`docs/barcs-input-output-examples.md`](docs/barcs-input-output-examples.md).
For the seven-method CRISPulator comparison and the real-data Waterbear
comparison, including the ranking-versus-calibration trade-off, use
[`docs/barcs-external-method-comparison.md`](docs/barcs-external-method-comparison.md).

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

- `R/bbreg.R`: dependency-free R implementation for one guide, contrasts,
  guide-by-guide screens, historical signed-score aggregation, exchangeable
  normal guide-beta inference, random-effects guide partial pooling, and
  empirical-Bayes heterogeneity moderation.
- `R/README.md`: source map, public workflow, and core input contracts.
- `examples/barcs_quickstart.R`: minimal longitudinal analysis with annotated
  counts, immutable library totals, sample metadata, and output.
- `julia/simulate_crispulator_facs.jl`: pinned CRISPulator 0.5.1 simulation of
  low 25%, high 25%, overlapping 0--100% bulk, and input samples.
- `examples/crispulator_facs_benchmark.R`: one-seed comparison of
  BARCS-original, BARCS-NORM, BARCS-partial, and BARCS-EB using one shared
  set of guide-level fits.
- `examples/crispulator_facs_repeated_benchmark.R`: the same four-method
  comparison over five seeds, MOI, guide quality, gene count, and replicate
  count.
- `examples/crispulator_facs_external_head_to_head.R`: one-evaluator
  comparison with official MAGeCK-MLE, edgeR-QL, DESeq2, and limma-voom,
  including result and input hashes.
- `examples/crispulator_facs_f1_threshold_curves.R`: standard gene-level F1,
  recall, precision, and realized-FDP curves at five nominal FDR thresholds.
- `docs/barcs-gene-methods.md`: equations, numerical interpretation,
  calibration, and diagnostics for the four guide-to-gene statistics.
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
  comparison of the same four BARCS gene statistics against the 26
  directionally validated genes and the selected 33-gene follow-up panel.
- `examples/waterbear_facs_external_head_to_head.R`: comparison of the four
  BARCS methods with rerun MAGeCK-MLE/MAGeCK test and the published
  Waterbear/MAUDE recovery totals.
- `examples/liang_cas13_benchmark.R`: five-cell-line Cas13 fitness
  processed-count sensitivity analysis. BARCS and MAGeCK-MLE fit a
  longitudinal slope across days 0, 7, and 14, with edgeR-QL, DESeq2, and
  limma-voom fitted to the same design.
- `examples/manuscript_liang_figure.R`: HAP1 guide-level volcano and
  longitudinal count-trajectory figure for the Liang analysis.
- `scripts/prepare_liang_cas13.R`,
  `scripts/count_liang_cas13_run.sh`, and
  `scripts/queue_liang_cas13_counts.sh`: download the Liang supplementary
  tables, stream the 30 longitudinal FASTQs through the published anchor/Bowtie
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
Rscript examples/barcs_quickstart.R
Rscript examples/cb2_generalization_proof.R
Rscript examples/cb2_generalization_walkthrough.R
Rscript examples/barcs_input_output_examples.R
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
julia --project=julia julia/simulate_crispulator_facs.jl
Rscript examples/crispulator_facs_repeated_benchmark.R
Rscript examples/crispulator_facs_external_head_to_head.R
Rscript examples/crispulator_facs_f1_threshold_curves.R
Rscript examples/waterbear_facs_benchmark.R
Rscript examples/waterbear_facs_external_head_to_head.R
Rscript examples/liang_cas13_benchmark.R
Rscript examples/manuscript_liang_figure.R
Rscript -e 'devtools::test("CB2")'
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

This branch intentionally uses only the CRISPulator FACS and Waterbear FACS
workflows to compare the four BARCS guide-to-gene methods. Other historical
benchmark scripts remain in the repository but do not define or tune these
four-method results.

The primary Liang comparison deliberately uses the deposited normalized guide
values so that the processed-data analysis is fully reproducible. The values
are fractional after median-ratio normalization, ComBat correction, and
outlier processing. The script rounds them once to the nearest pseudo-count
and supplies the identical matrix to all newly fitted methods. Accordingly,
this is a processed-count sensitivity analysis: none of the count-based
methods retains its literal raw-count sampling interpretation.
Known-essential protein-coding genes are positives, and cell-line-specific
non-expressed lncRNAs are null controls.

All compared methods use the three deposited time points and estimate a
continuous time slope. A replicate block is included for complete
trajectories; K562 lacks one baseline replicate and therefore uses an
unblocked slope. The longitudinal regression comparisons use two-sided
significance values and FDRs.

Adding day 7 raises BARCS macro-average precision from 0.776 to 0.786 and
FDR-0.10 essential-gene recall from 0.487 to 0.580. BARCS has lower null
calibration error than longitudinal MAGeCK-MLE (0.035 versus 0.103), although
MAGeCK-MLE has higher average precision (0.833) and FDR-thresholded recall
(0.690). This is the intended boundary of the claim: BARCS generalizes
beta-binomial inference to longitudinal regression and improves calibration
relative to the matched negative-binomial regression without claiming
an endpoint-only analysis.

BARCS-original, edgeR-QL, DESeq2, and limma-voom all use unweighted signed-\(z\)
aggregation of two-sided guide probabilities after their respective
guide-level fits. This holds the gene combiner fixed while comparing
guide-level models. MAGeCK-MLE instead uses its native joint gene model.

The primary Liang result uses every valid guide. Mean BARCS--MAGeCK-MLE
effect-rank correlation is 0.875 across cell lines. Versioned metrics,
concordance estimates, and the rounding audit are under `data/derived/`.

An optional raw-read confirmation remains available. It streams all
longitudinal runs without retaining FASTQs:

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

The four-method benchmark does not require MAGeCK, edgeR, DESeq2, or limma.
All four methods reuse the same dependency-free BARCS guide fits. The
RcppArmadillo weighted-crossproduct kernels from the pinned CB2 submodule are
used automatically when the compiled library is available.

The GSE70038 script expects official MAGeCK 0.5.9.5 at
`.venv/bin/mageck`. A project-local installation can be prepared with:

```sh
python3 -m venv .venv
.venv/bin/pip install numpy scipy
curl -L -o /tmp/mageck-0.5.9.5.tar.gz \
  'https://sourceforge.net/projects/mageck/files/0.5/mageck-0.5.9.5.tar.gz/download'
.venv/bin/pip install /tmp/mageck-0.5.9.5.tar.gz
```

The four-method CRISPulator analysis contains ten one-at-a-time parameter
scenarios, five fixed seeds per scenario, and one shared guide-level fit for
all four gene statistics. Across the 50 low--bulk--high runs,
BARCS-original has the strongest average ranking and threshold power:
average precision 0.856, AUROC 0.918, directional recall 0.638, and F1 0.702.
Its mean realized FDP is 0.085 and 5.5% of negative-control gene p-values fall
below 0.05.

BARCS-EB trades power for calibration. Its average precision is 0.840,
directional recall 0.463, and F1 0.599, while realized FDP falls to 0.013 and
the negative-control p-value frequency falls to 0.013. BARCS-partial is
weaker overall: average precision 0.779, directional recall 0.390, F1 0.525,
and realized FDP 0.052. BARCS-original has the highest mean average precision
in nine of ten scenarios; BARCS-EB leads one. Therefore the current results do
not support replacing the historical statistic universally.

BARCS-NORM implements the exchangeable normal guide-beta model directly. Its
average precision is 0.706, directional recall is 0.159, F1 is 0.246, and
realized FDP is 0.059. Its negative-control \(p<0.05\) rate is well calibrated
at 0.047, but estimating a separate standard deviation from roughly five
guides leaves only about four reference degrees of freedom and sharply limits
power. A plug-in standard-normal reference is also reported by the function
but is not used for primary calls because it treats the estimated standard
deviation as known.

Against external methods across the same 50 simulations, edgeR-QL has the
highest mean average precision (0.877) and directional recall (0.822), while
limma-voom has the highest F1 (0.772). Their nominal FDR 0.10 thresholds are
anti-conservative in this simulation: realized FDP is 0.277 for edgeR-QL,
0.233 for DESeq2, and 0.257 for limma-voom. BARCS-original has lower average
precision (0.856) and F1 (0.702), but a realized FDP of 0.085 and the
negative-control p-value rate closest to 0.05 (0.055). MAGeCK-MLE is the
closest external compromise, with average precision 0.849, F1 0.683, and
realized FDP 0.058. The complete interpretation and provenance audit are in
`docs/barcs-external-method-comparison.md`.

The FDR-threshold scan shows that this is not only a nominal-threshold
artifact. Across all 50 runs, edgeR-QL at nominal FDR 0.01 has mean F1 0.796
and realized FDP 0.087, compared with BARCS-original at nominal FDR 0.10 with
F1 0.702 and realized FDP 0.085. This post hoc matched-FDP observation
identifies a real guide-level information-borrowing advantage in the
simulation; it does not establish 0.01 as a prospectively calibrated edgeR
threshold for other screens.

At the four-replicate baseline, mean average precision is 0.917, 0.777, 0.842,
and 0.903 for original, NORM, partial, and EB, respectively. Their directional
recalls are 0.764, 0.209, 0.450, and 0.551, while realized FDPs are 0.086,
0.049, 0.046, and 0.008.
In the diagnostic one-replicate setting, original makes no FDR 0.10 calls;
partial and EB reach directional recalls 0.128 and 0.117, respectively.
These remain reagent-consistency results rather than biological-replicate
inference.

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

The GSE242880 comparison uses the low-coverage, high-MOI primary-T-cell arm:
four ordered IL2RA FACS bins for each of three donors. The 26 evaluation genes
were validated by individual knockout and flow cytometry. BARCS-original,
BARCS-NORM, BARCS-partial, and BARCS-EB recover 22, 0, 19, and 23 of 26 genes
in the expected direction at gene FDR 0.10, with 49, 0, 71, and 60 total
screen discoveries.

The deposited follow-up table also contains seven candidates that did not
validate experimentally. In this selected 33-gene panel, original, partial,
and EB have F1 values 0.863, 0.792, and 0.885 and balanced accuracies 0.709,
0.651, and 0.728; NORM has F1 0 because it makes no calls. Original retains
the highest average precision (0.945
versus 0.940 for EB), whereas EB has the highest validated recovery, F1, and
balanced accuracy. The selected panel is supporting evidence, not an unbiased
genome-wide negative set.

BARCS-NORM still ranks the selected panel above chance (average precision
0.897) but makes no FDR-0.10 calls, again showing the small-guide
standard-deviation penalty. In the external GSE242880 comparison, MAGeCK-MLE
recovers 17/26 validated
genes with 72 calls and outer-bin MAGeCK test recovers 18/26 with 30 calls.
Published Waterbear and MAUDE totals are 24/26 with 79 calls and 25/26 with
406 calls. Those published aggregates do not provide complete per-gene
scores for the selected 33-gene panel, so F1 and average precision are
reported only for the six methods rerun from complete outputs.

The shared raw guide fit is inflated among the 593 non-targeting guides:
13.3% have nominal p-values below 0.05. Historical control calibration reduces
that fraction to 5.1%. Partial and EB instead use gene-statistic calibration;
because all non-targeting guides share one deposited gene label, their
gene-level null falls back to the robust whole-screen center and scale.
FACS bins remain correlated partitions, so none of these methods replaces a
specialist joint-bin model. Compact results are versioned under
the legacy `data/derived/waterbear_facs_three_method_*` filenames.

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

Within BARCS, the current evidence supports BARCS-original when ranking and
recall are primary, and BARCS-EB when conservative null calibration is
primary. BARCS-NORM is the most literal unweighted normal model of guide beta
values, but its per-gene variance estimate is underpowered with three to six
guides. BARCS-partial exposes guide heterogeneity and influence diagnostics
but is not the best default in the present benchmarks. This choice must be
made from the inferential goal or an external protocol, not selected
retrospectively from whichever method gives the most favorable result.

BARCS is not a model for a continuous phenotype measured per cell when guide
identity and phenotype are not jointly observed. The current regression layer
is a transparent research implementation. The guide hierarchy does not turn
multiple reagents into biological replicates, and repeated bins or donors
still require an appropriate dependence model.
