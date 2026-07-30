# BARCS

**B**eta-binomial **A**nalysis and **R**egression for **C**RISPR **S**creens.

BARCS fits a beta-binomial regression to pooled CRISPR screen counts on an
arbitrary design matrix, so a screen can be analyzed as a *coefficient*
question rather than a two-group comparison: does guide abundance follow dose,
time, an ordered phenotype, a donor-adjusted treatment effect, or a
knockout-by-treatment interaction?

This repository holds three things:

- **the method** — a dependency-free R implementation in [`R/bbreg.R`](R/bbreg.R),
  packaged as `BARCS`;
- **the manuscript** — [`main.tex`](main.tex) and [`sections/`](sections/), the
  derivation, benchmarks, and limitations;
- **the benchmark suite** — [`examples/`](examples/) and
  [`scripts/`](scripts/), the scripts that produced every number in the
  manuscript, with their outputs versioned under
  [`data/derived/`](data/derived/).

## Why a regression, and why beta-binomial

CB² ([Jeong et al. 2019](https://doi.org/10.1101/gr.245571.118)) answers a
two-condition question: does guide abundance differ between two groups? That
is the wrong shape for a longitudinal, donor-blocked, or factorial screen.
BARCS keeps CB²'s sampling idea and replaces the two-group mean with an
ordinary R model matrix.

A guide's count $Y_i$ in library $i$ is modeled conditional on the
**unfiltered** full-library total $N_i$. Two sources of variation are kept
separate:

$$\operatorname{Var}\!\left(\frac{Y_i}{N_i}\right)
  = \mu(1-\mu)\left\{\frac{1-\rho}{N_i} + \rho\right\}$$

The first term is sequencing precision, which shrinks as you sequence deeper.
The second is between-library heterogeneity, which does not. The intraclass
correlation $\rho$ is estimated per guide from a Pearson moment equation, the
mean follows a logit link on the design matrix, and coefficients are fitted by
feasible iteratively reweighted least squares with weights
$w_i = N_i\mu(1-\mu)/\{1 + (N_i-1)\rho\}$.

The practical consequence: a coefficient is referred to $t_{m-q}$ for $m$
independently sequenced libraries and $q$ fitted parameters — **degrees of
freedom come from libraries, not from reads.** Sequencing a screen more deeply
does not buy inferential confidence about a sample-level effect.

Three data contracts follow from this and are enforced by the code:

| Object | Meaning | Common mistake |
|---|---|---|
| `counts` | guide × library count matrix | — |
| `totals` | unfiltered mapped-guide total per library | recomputing it *after* filtering guides, which shrinks the denominator and inflates significance |
| `data` | one row per library, in column order | passing per-cell or per-replicate metadata instead |

## Quickstart

The reference implementation is base R only, so it runs without installing
anything:

```r
source("R/bbreg.R")

counts <- rbind(
  guide_1 = c(120, 110, 85, 70, 48, 33),
  guide_2 = c( 75,  79, 72, 70, 68, 66),
  control = c( 42,  45, 41, 43, 44, 40)
)
colnames(counts) <- paste0("sample_", 1:6)

# All mapped guides per library — not the column sums of the filtered matrix.
library_totals <- c(100000, 98000, 105000, 101000, 99000, 103000)

sample_data <- data.frame(
  sample    = colnames(counts),
  replicate = factor(rep(c("A", "B"), times = 3)),
  day       = rep(c(0, 7, 14), each = 2)
)

result <- bb_screen(
  counts  = counts,
  totals  = library_totals,
  data    = sample_data,
  formula = ~ replicate + I(day / 14),
  term    = "I(day/14)",
  guide   = rownames(counts)
)
```

Runnable version: [`examples/barcs_quickstart.R`](examples/barcs_quickstart.R).
Printed input → call → output tables:
[`docs/barcs-input-output-examples.md`](docs/barcs-input-output-examples.md).

### The public API

| Function | Purpose |
|---|---|
| `bbreg()` | fit one guide on a design matrix |
| `bb_contrast()` | test a named linear combination of coefficients |
| `bb_screen()` | apply the fit across a count matrix; one tidy row per guide |
| `bb_calibrate_controls()` | estimate an empirical null scale from prespecified negative controls |
| `bb_moderate_dispersion()` | limma-style empirical-Bayes moderation of guide dispersion |
| `bb_gene_*()` | optional, explicitly labeled guide-to-gene summaries |

Source map and internals: [`R/README.md`](R/README.md).

To install as a package (needs the `CB2` submodule, see below):

```sh
R CMD INSTALL CB2 && R CMD INSTALL .
```

The compiled RcppArmadillo kernels in [`src/`](src/) are used automatically when
available; the R path is the reference and remains correct without them.

## Relationship to CB²

BARCS is an additive path, not a replacement for the two-group workflow.

| | CB² | BARCS |
|---|---|---|
| Primary question | difference between two groups | trend, adjusted association, interaction, or contrast |
| Design | reference/comparison labels | arbitrary full-rank model matrix |
| Effect | group difference or fold change | conditional log-odds coefficient |
| Dispersion | estimated within groups | one guide-wise $\rho$ across the fitted design |
| Estimator | weighted group proportions | feasible logistic IRLS |
| Reference df | Welch–Satterthwaite | residual sample df |
| Variance principle | beta-binomial | same beta-binomial principle |

The two share a library-total-conditional variance principle and nothing more.
An identity-link model can reproduce already-computed CB² group summaries when
its weights are defined from those summaries, but that is a software check, not
a nesting theorem. **BARCS makes no formal equivalence claim.**

## What the benchmarks actually show

The honest summary is that BARCS makes a class of designs analyzable under a
beta-binomial sampling model. It does **not** win on ranking or calibration
against mature negative-binomial and RNA-seq methods. Every number below is
reproduced by a script in `examples/` and deposited under `data/derived/` or
`results/`.

| Benchmark | Result |
|---|---|
| **Liang Cas13**, 4 cell lines, longitudinal slope | After applying the *same* aggregation-matched control-scaling rule to all five methods, mean absolute calibration error was 0.0171 (MAGeCK-MLE), 0.0178 (limma-voom), 0.0181 (edgeR-QL), 0.0198 (DESeq2), 0.0207 (BARCS). BARCS was last on average precision (0.838 vs 0.874–0.877) and on FDR-0.10 recall (0.600 vs 0.738–0.767). Control scaling is method-agnostic, not a beta-binomial advantage. |
| **GSE242880 IL2RA**, ordered FACS bins, donor-adjusted | 22/26 validated regulators at 49 discoveries (44.9% yield); held-out non-targeting rate 0.049 after five-fold cross-fitting, versus 0.133 raw. The outer-bin fit was *more* efficient per call (52.9%), and Waterbear (24/26, 79 calls) and MAUDE (25/26, 406 calls) recovered more. Extra bins buy sensitivity, not precision. |
| **HT-29 serial harvests** (Tzelepis) | AUROC 0.9786 (BARCS time), 0.9789 (MAGeCK time), 0.9747 (Chronos) — and 0.9792 for the deposited *day-25 endpoint*. The trajectory offers no measured ranking advantage. Reported as descriptive: seven harvests of one lineage are not seven biological replicates. |
| **Sanson A375** endpoint | AUROC 0.9598 for both BARCS and MAGeCK-MLE; after CNV correction the difference is 0.00316 (95% CI −0.00265 to 0.00924). No advantage for either count model. |
| **GSE70038**, four terminal coefficients jointly | Top-200 depletion lists share 139–144 genes; exclusive hits are essential-process-coherent for both methods (62.3–75.9% BARCS-only, 59.0–86.2% MAGeCK-only against the Hart reference). Complementary, not a ranking. |
| **CRISPulator**, 10,000 genes, MOI 0.20 | Dispersion moderation helps within BARCS: average precision 0.902 → 0.919, F1 0.811 → 0.847, realized FDP unchanged (0.065 → 0.066). MAGeCK-MLE reached 0.921 AP at FDP 0.006; CRISPhieRmix 0.867 at FDP 0.197. Fitting 50,000 guide regressions took 70.1–102.8 s per seed. |
| **simCRISPR** knockout × treatment interaction | A full-library denominator shifted true-zero guides (median fitted interaction 0.184) because widespread depletion raises every survivor's library share. A control denominator removed the shift and improved F1 in *every* seed — approximately 0.71 → 0.91 in the two runs that did not collapse. |
| **External false-discovery audit** (Dempster et al.) | The reported 152 CB² discoveries among Avana null genes fell to 0 once full-library totals were restored, i.e. the count was a consequence of pre-fit guide subsetting. Nominal-level calibration remains unresolved: the largest cell-line-specific null $p<0.05$ rate was still 0.145. Disclosed as an audit of the author's own prior method. |

## Where BARCS is weakest

Read this before using it for a confirmatory analysis. These are measured
limitations, not hypothetical ones.

**Gene-level rates can be anti-conservative when guides are correlated.** The
guide-to-gene step combines guide $p$-values with a directional Stouffer rule
under a $\sqrt{m_g}$ independence reference. In an all-null grid
(`examples/barcs_null_calibration_grid.R`), imposing within-gene guide
correlation $r = 0.4$ left *guide*-level type-I error at or below 0.048 while
raising *gene*-level error to 0.108–0.144 (3 guides/gene) and 0.224–0.262
(5 guides/gene). Aggregation-matched split-control scaling reduced this to
0.056–0.112 but did not remove it. Real guides against one gene share target
biology, efficiency, and seed-based off-target effects, so treat published
gene-level FDRs as possibly ~2× anti-conservative for five-guide genes. A
hierarchical gene model is the right fix and is not implemented.

**The covariance is a plug-in.** $(X'WX)^{-1}$ treats $\hat\rho$ as known and
the $t$ reference does not propagate dispersion-estimation uncertainty, so
small-sample calibration is not guaranteed. In the committed continuous-dose
simulation, unmoderated BARCS had gene-level type-I error 0.081 and realized
FDP 0.130 at nominal 0.05; moderation cut FDP to 0.091 without repairing the
marginal error. Official MAGeCK-MLE was 0.038 / 0.049 in the same realization.
Expanding dispersion over the full range did *not* reproduce this failure —
the gene combiner above is the larger contributor.

**Control calibration must be held out.** `bb_calibrate_controls()` estimates
its scale from the control tail, so evaluating the same guides is circular. Use
cross-fitting or a disjoint control split, and estimate the scale at the same
aggregation level as the statistic you report. Controls must also share one
design.

**$\hat\rho$ is truncated at zero** when the Pearson statistic falls below its
degrees of freedom. Observed frequency in real fits: 11/86,840 (HT-29),
0/5,999 (IL2RA), 1,738/280,541 (Liang) —
see `examples/barcs_real_data_boundary_audit.R`.

**Guides are not biological replicates**, and correlated partitions of one cell
pool are not independent libraries. For several FACS bins from one pool, use a
specialist joint model such as Waterbear; BARCS's negative controls diagnose the
independence violation and calibrate one operating point, but they do not make
the bin margins independent. BARCS is also not a model for a per-cell
continuous phenotype when guide identity and phenotype are not jointly
observed.

## Repository layout

The R package is developed in
[`jeonglab-bcm/CB2`](https://github.com/jeonglab-bcm/CB2) and tracked here as
the `CB2/` submodule, so the parent repository records the exact package commit
behind every manuscript revision.

```sh
git clone --recurse-submodules https://github.com/jeonglab-bcm/BARCS.git
cd BARCS
```

Package changes are committed from inside `CB2/`; the parent then commits the
updated submodule pointer. Full two-repository workflow:
[`DEVELOPMENT.md`](DEVELOPMENT.md).

| Path | Contents |
|---|---|
| `R/bbreg.R` | the implementation: single-guide fit, contrasts, screens, control calibration, dispersion moderation, guide-to-gene summaries |
| `src/` | RcppArmadillo weighted-crossproduct and symmetric-solve kernels |
| `main.tex`, `sections/` | manuscript source |
| `output/pdf/` | rendered manuscript |
| `examples/` | every benchmark and simulation script |
| `julia/` | pinned CRISPulator 0.5.1 FACS simulation |
| `scripts/` | Liang data preparation, MAGeCK 0.5.9.5 compatibility and CNV shims |
| `data/derived/` | versioned benchmark tables, per-seed values, and audits |
| `results/` | generated benchmark intermediates (gitignored except deposited summaries) |
| `docs/` | gene-summary methods, input/output examples, external comparison, peer-review records |
| `tests/run_tests.R` | base-R regression and input-validation tests |

Selected scripts:

- `examples/liang_cas13_benchmark.R` — four-cell-line Cas13 longitudinal slope
  against MAGeCK-MLE, edgeR-QL, DESeq2, limma-voom on one shared matrix
- `examples/waterbear_facs_benchmark.R` — ordered four-bin IL2RA analysis with
  cross-fitted control calibration
- `examples/crispulator_facs_moi_10k_benchmark.R` — genome-scale simulation with
  known truth
- `examples/simcrispr_interaction_benchmark.R` — factorial interaction and
  denominator ablation
- `examples/barcs_null_calibration_grid.R` — the all-null calibration grid,
  including the correlated-guide arm
- `examples/simulation.R` — the prespecified continuous-dose diagnostic
- `examples/sanson_benchmark.R`, `examples/chronos_tzelepis_benchmark.R`,
  `examples/gse70038_comparison.R` — real-data endpoint, time-course, and
  multi-coefficient comparisons

## Reproducing

```sh
Rscript tests/run_tests.R
Rscript examples/barcs_quickstart.R
Rscript examples/simulation.R
Rscript examples/barcs_null_calibration_grid.R

# Genome-scale simulation (Julia + CRISPulator, pinned)
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
julia --project=julia julia/simulate_crispulator_facs.jl
Rscript examples/crispulator_facs_moi_10k_benchmark.R

# Real-data benchmarks
Rscript examples/liang_cas13_benchmark.R
Rscript examples/waterbear_facs_benchmark.R
Rscript examples/sanson_benchmark.R

# Manuscript
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

External comparators are kept external on purpose — no partial
reimplementations. MAGeCK comparisons expect official MAGeCK 0.5.9.5 at
`.venv/bin/mageck`:

```sh
python3 -m venv .venv
.venv/bin/pip install numpy scipy
curl -L -o /tmp/mageck-0.5.9.5.tar.gz \
  'https://sourceforge.net/projects/mageck/files/0.5/mageck-0.5.9.5.tar.gz/download'
.venv/bin/pip install /tmp/mageck-0.5.9.5.tar.gz
```

The Liang analysis deliberately uses the deposited normalized values so the
processed-data comparison is fully reproducible. Those values are fractional
after median-ratio normalization, ComBat correction, and outlier processing; the
script rounds them once and gives the identical matrix to every method. It is
therefore a processed-count sensitivity analysis, and no count-based method
retains its literal raw-count sampling interpretation. An optional raw-read
confirmation streams the FASTQs without retaining them:

```sh
bash scripts/queue_liang_cas13_counts.sh 2
pueue wait --group liang-cas13
```

## Scope

Use CB² for a direct two-condition comparison. Use a specialist joint model
when several bins are correlated partitions of one biological pool. Use BARCS
when each independently sequenced library carries a quantitative or
multivariable sample-level design — and choose the guide-to-gene summary from
the inferential goal or a preregistered protocol, never retrospectively from
whichever one gives the most favorable result.

BARCS is a transparent research implementation. Confirmatory use requires
independent biological replication, likelihood-compatible counts, and held-out
or cross-fitted controls.

## Citation

The two-group special case and its method:

> Jeong H-H, Kim SY, Rousseaux MWC, Zoghbi HY, Liu Z. Beta-binomial modeling of
> CRISPR pooled screen data identifies target genes with greater sensitivity and
> fewer false negatives. *Genome Research* 29:999–1008 (2019).
> https://doi.org/10.1101/gr.245571.118

The overdispersed regression construction follows Williams (1982) and
[Baggerly et al. (2004)](https://doi.org/10.1186/1471-2105-5-144). The BARCS
manuscript is in [`main.tex`](main.tex) and rendered under
[`output/pdf/`](output/pdf/).

## License

MIT — see [LICENSE](LICENSE).
