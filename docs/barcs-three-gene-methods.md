# Three BARCS guide-to-gene methods

All three methods begin with the same guide-by-sample count matrix and the
same beta-binomial regression. They differ only after one phenotype
coefficient and standard error have been estimated for every guide.

For a gene with four guides, the shared input might be:

| guide | phenotype coefficient | standard error | two-sided p-value |
|---|---:|---:|---:|
| g1 | -0.75 | 0.20 | 0.00018 |
| g2 | -0.90 | 0.25 | 0.00032 |
| g3 | -0.60 | 0.30 | 0.04550 |
| g4 | 0.10 | 0.40 | 0.80259 |

No method treats these four guides as four biological samples.

## 1. BARCS-original

The historical calculation converts every guide p-value to a signed normal
score,

\[
z_{gj}
=
\operatorname{sign}(\widehat b_{gj})
\Phi^{-1}(1-p_{gj}/2),
\]

and calculates

\[
Z_g=\frac{\sum_j z_{gj}}{\sqrt{m_g}}.
\]

This method uses each guide's direction and significance, but it does not
estimate an explicit guide-disagreement variance.

## 2. BARCS-partial

The partial-pooling model works with guide coefficients rather than guide
p-values:

\[
\widehat b_{gj}\mid\beta_g,\tau_g^2
\sim
N\left(\beta_g,\,
\operatorname{SE}(\widehat b_{gj})^2+\tau_g^2\right).
\]

Here, \(\beta_g\) is the shared gene effect and \(\tau_g^2\) measures how much
the guide effects disagree beyond their estimated sampling errors.

The raw disagreement estimate is the non-negative
DerSimonian--Laird estimate,

\[
\widehat\tau_g^2
=
\max\left\{0,\frac{Q_g-(m_g-1)}{C_g}\right\}.
\]

The gene effect is then

\[
\widehat\beta_g
=
\frac{\sum_j w_{gj}\widehat b_{gj}}{\sum_jw_{gj}},
\qquad
w_{gj}
=
\frac{1}
{\operatorname{SE}(\widehat b_{gj})^2+\widehat\tau_g^2}.
\]

Consistent guides give \(\widehat\tau_g^2\) close to zero. Contradictory guides
increase \(\widehat\tau_g^2\), reduce their weights, and increase the gene
standard error.

## 3. BARCS-EB

With only three to six guides per gene, each
\(\widehat\tau_g^2\) is noisy. BARCS-EB first estimates the typical
screen-wide guide-disagreement variance, \(\tau_0^2\), from pooled excess
Cochran \(Q\). It then moderates each gene:

\[
\widetilde\tau_g^2
=
\frac{
(m_g-1)\widehat\tau_g^2+d_0\tau_0^2
}{
(m_g-1)+d_0
}.
\]

The current prespecified value is \(d_0=4\). A typical five-guide gene also
has four heterogeneity degrees of freedom, so its raw estimate and the
screen-wide estimate receive equal nominal weight. The prior scale
\(\tau_0^2\) is estimated from the screen without using simulated truth or
validation labels.

BARCS-EB uses \(\widetilde\tau_g^2\) in the same weighted gene-effect formula
as partial pooling. It is therefore a stabilized version of BARCS-partial,
not a fourth guide-level regression.

## Calibration and diagnostics

BARCS-original reproduces the historical guide-level control calibration
followed by signed-score aggregation. BARCS-partial and BARCS-EB calibrate
their gene statistics against negative-control genes when at least ten are
available and otherwise use the robust center and scale of the whole
gene-statistic distribution.

The two guide-effect methods report:

- raw and moderated guide-disagreement variance;
- \(I^2\), the fraction of observed guide variation attributed to
  disagreement;
- guide-direction agreement;
- the largest guide weight;
- the largest leave-one-guide-out change;
- whether the gene-effect direction survives every leave-one-guide-out fit.

These diagnostics preserve information that would disappear if guide counts
were summed before fitting.

## Benchmark scope

The development branch evaluates these three gene statistics only in:

1. `examples/crispulator_facs_benchmark.R`;
2. `examples/crispulator_facs_repeated_benchmark.R`, which creates the
   CRISPulator multi-parameter, five-seed, and replicate summaries;
3. `examples/waterbear_facs_benchmark.R`.

Every within-benchmark comparison reuses the identical guide-level fit.
