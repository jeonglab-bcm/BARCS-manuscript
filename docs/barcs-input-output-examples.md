# BARCS input and output examples

BARCS stands for Beta-binomial Analysis and Regression for CRISPR Screens.

Run every example from the repository root with:

```sh
Rscript examples/barcs_input_output_examples.R
```

The runnable source is
[`examples/barcs_input_output_examples.R`](../examples/barcs_input_output_examples.R).

The original-CB² calculation below is a legacy compatibility example. It
shows that saturated GLS reproduces completed CB² group summaries; it does not
claim that default `bbreg(~ group)` is the same finite-sample estimator.

## Minimum input structure

For one guide, `bbreg()` needs:

- `count`: the guide read count in each independently sequenced sample;
- `total`: the total mapped guide reads in the corresponding sample;
- `data`: one row of sample information per count;
- `formula`: the columns of `data` to include in the model.

The ordering must match:

$$
(K_1,\ldots,K_m)
\longleftrightarrow
(N_1,\ldots,N_m)
\longleftrightarrow
\text{rows }1,\ldots,m\text{ of the sample table}.
$$

## Example 1: two groups

### Input

| Sample | Group | Guide count | Total reads | Observed proportion |
|---|---|---:|---:|---:|
| A1 | A | 74 | 52,000 | 0.001423077 |
| A2 | A | 112 | 81,000 | 0.001382716 |
| A3 | A | 91 | 69,000 | 0.001318841 |
| A4 | A | 139 | 97,000 | 0.001432990 |
| B1 | B | 128 | 55,000 | 0.002327273 |
| B2 | B | 177 | 76,000 | 0.002328947 |
| B3 | B | 164 | 72,000 | 0.002277778 |
| B4 | B | 211 | 93,000 | 0.002268817 |

```r
two_group <- data.frame(
  group = factor(rep(c("A", "B"), each = 4)),
  count = c(74, 112, 91, 139, 128, 177, 164, 211),
  total = c(52000, 81000, 69000, 97000, 55000, 76000, 72000, 93000)
)
```

### Input to original CB²

Original CB² fits each group separately:

```r
fit_a <- fit_ab(
  matrix(two_group$count[two_group$group == "A"], nrow = 1),
  matrix(two_group$total[two_group$group == "A"], nrow = 1)
)

fit_b <- fit_ab(
  matrix(two_group$count[two_group$group == "B"], nrow = 1),
  matrix(two_group$total[two_group$group == "B"], nrow = 1)
)
```

It tests

$$
T_{\mathrm{CB2}}
=
\frac{\widehat p_B-\widehat p_A}{\sqrt{V_A+V_B}}.
$$

### Original CB² output

| $\widehat p_A$ | $\widehat p_B$ | Effect $B-A$ | SE | t | df | p |
|---:|---:|---:|---:|---:|---:|---:|
| 0.001391304 | 0.002297297 | 0.000905993 | 0.000111311 | 8.139330 | 5.647255 | 0.000251293 |

### Input to default BARCS

```r
group_fit <- bbreg(
  count = two_group$count,
  total = two_group$total,
  formula = ~ group,
  data = two_group
)
```

The model is

$$
\operatorname{logit}(\mu_i)
=\beta_0+\beta_{\mathrm{group}}I(\mathrm{group}_i=B).
$$

### Default BARCS output

| Coefficient | Estimate | SE | t | df | p |
|---|---:|---:|---:|---:|---:|
| Intercept | -6.576121319 | 0.049063177 | -134.033745 | 6 | $1.16\times10^{-11}$ |
| Group B | 0.502399324 | 0.062299034 | 8.064320 | 6 | 0.000194616 |

The original effect is a raw-proportion difference. The BARCS effect is a
log-odds difference:

$$
\widehat\beta_{\mathrm{group}}
=\operatorname{logit}(\widehat\mu_B)
-\operatorname{logit}(\widehat\mu_A).
$$

Consequently, the effect numbers have different units. Their t-statistics are
close, but finite-sample identity is not claimed.

The model-based BARCS standard errors shown below treat the guide-wise
dispersion estimate as a fixed plug-in value. Residual Student t degrees of
freedom do not formally account for dispersion-estimation uncertainty.

## Example 2: continuous dose

### Input

| Dose | Total reads | Guide count | Observed proportion |
|---:|---:|---:|---:|
| 0 | 80,000 | 73 | 0.0009125 |
| 0 | 100,000 | 91 | 0.0009100 |
| 1 | 80,000 | 103 | 0.0012875 |
| 1 | 100,000 | 129 | 0.0012900 |
| 2 | 80,000 | 147 | 0.0018375 |
| 2 | 100,000 | 183 | 0.0018300 |
| 3 | 80,000 | 208 | 0.0026000 |
| 3 | 100,000 | 260 | 0.0026000 |
| 4 | 80,000 | 295 | 0.0036875 |
| 4 | 100,000 | 368 | 0.0036800 |

```r
dose_fit <- bbreg(
  count = dose_input$count,
  total = dose_input$total,
  formula = ~ dose,
  data = dose_input
)
```

The model is

$$
\operatorname{logit}(\mu_i)
=\beta_0+\beta_{\mathrm{dose}}d_i.
$$

### Output

| Coefficient | Estimate | SE | t | df | p |
|---|---:|---:|---:|---:|---:|
| Intercept | -7.000957876 | 0.052622887 | -133.040171 | 8 | $1.14\times10^{-14}$ |
| Dose | 0.350250934 | 0.017728616 | 19.756248 | 8 | $4.49\times10^{-8}$ |

The data were created with a true dose coefficient of $0.35$. The fitted value
is

$$
\widehat\beta_{\mathrm{dose}}=0.350251.
$$

Its odds-ratio interpretation is

$$
\exp(\widehat\beta_{\mathrm{dose}})
=\exp(0.350251)
\approx1.419.
$$

Thus each one-step dose increase is associated with approximately a 41.9%
increase in the guide-abundance odds.

## Example 3: continuous dose with batch adjustment

### Input

| Sample | Dose | Batch | Total reads | Guide count | Observed proportion |
|---|---:|---|---:|---:|---:|
| s01 | 0 | A | 80,000 | 73 | 0.0009125 |
| s02 | 1 | A | 95,000 | 117 | 0.0012316 |
| s03 | 2 | A | 85,000 | 141 | 0.0016588 |
| s04 | 3 | A | 105,000 | 235 | 0.0022381 |
| s05 | 4 | A | 90,000 | 272 | 0.0030222 |
| s06 | 5 | A | 100,000 | 407 | 0.0040700 |
| s07 | 0 | B | 80,000 | 89 | 0.0011125 |
| s08 | 1 | B | 95,000 | 143 | 0.0015053 |
| s09 | 2 | B | 85,000 | 172 | 0.0020235 |
| s10 | 3 | B | 105,000 | 287 | 0.0027333 |
| s11 | 4 | B | 90,000 | 332 | 0.0036889 |
| s12 | 5 | B | 100,000 | 497 | 0.0049700 |

```r
adjusted_fit <- bbreg(
  count = adjusted_input$count,
  total = adjusted_input$total,
  formula = ~ dose + batch,
  data = adjusted_input
)
```

The model is

$$
\operatorname{logit}(\mu_i)
=\beta_0+\beta_{\mathrm{dose}}d_i
+\beta_{\mathrm{batchB}}I(\mathrm{batch}_i=B).
$$

### Output

| Coefficient | Estimate | SE | t | df | p |
|---|---:|---:|---:|---:|---:|
| Intercept | -6.999085429 | 0.050551579 | -138.454339 | 9 | $2.72\times10^{-16}$ |
| Dose | 0.299891857 | 0.012326192 | 24.329644 | 9 | $1.60\times10^{-9}$ |
| Batch B | 0.200188653 | 0.038282790 | 5.229207 | 9 | 0.000542500 |

The data were created with

$$
\beta_{\mathrm{dose}}=0.3,
\qquad
\beta_{\mathrm{batchB}}=0.2.
$$

The fitted values recover both while keeping the dose effect adjusted for
batch.

### Output for a two-dose-step contrast

```r
bb_contrast(adjusted_fit, c(dose = 2))
```

This tests

$$
2\beta_{\mathrm{dose}}=0.
$$

| Estimate | SE | t | df | p |
|---:|---:|---:|---:|---:|
| 0.599783714 | 0.024652383 | 24.329644 | 9 | $1.60\times10^{-9}$ |

Doubling a coefficient and its standard error does not change its t-statistic:

$$
\frac{2\widehat\beta_{\mathrm{dose}}}
{2\operatorname{SE}(\widehat\beta_{\mathrm{dose}})}
=
\frac{\widehat\beta_{\mathrm{dose}}}
{\operatorname{SE}(\widehat\beta_{\mathrm{dose}})}.
$$

## Example 4: guide consistency when the screen has one replicate

Suppose an irreplaceable screen has only one low, bulk, and high library, but
five independently designed guides target each gene. The sample-level design
still has only one residual degree of freedom. The guides therefore cannot be
called biological replicates.

`bb_gene_consistency()` instead asks a narrower question: do the guide
coefficient estimates support a shared gene effect more strongly than the
empirical null?
Its input is a guide-level result with at least `gene`, `estimate`, and
`std_error` columns:

```r
guide_consistency <- bb_gene_consistency(
  guide_results,
  control = grepl("^NTC", guide_results$gene)
)
```

For guide \(j\) targeting gene \(g\), let
\(w_{gj}=\operatorname{SE}(\widehat\beta_{gj})^{-2}\). The function directly
estimates a shared gene coefficient and its model-based standard error:

$$
\widehat\beta_g
=
\frac{\sum_j w_{gj}\widehat\beta_{gj}}{\sum_jw_{gj}},
\qquad
\operatorname{SE}(\widehat\beta_g)
=
\frac{1}{\sqrt{\sum_jw_{gj}}},
\qquad
T_g
=
\frac{\widehat\beta_g}{\operatorname{SE}(\widehat\beta_g)}.
$$

No guide p-values are combined by Fisher or Stouffer.
Ten five-guide control genes set the null center and contribute a tail-scale
check; the all-gene median absolute deviation supplies a robust scale. The
larger scale is used, with a lower bound of one. For the deterministic example,
the null center is approximately zero and the scale is one.

| Gene | Guides | Shared effect | SE | Raw Wald score | Calibrated z | Direction agreement | p | FDR |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| null_gene | 5 | -0.002 | 0.0447 | -0.0447 | -0.0447 | 0.50 | 0.964 | 0.964 |
| signal_gene | 5 | 0.584 | 0.0447 | 13.0586 | 13.0586 | 1.00 | \(5.67\times10^{-39}\) | \(6.81\times10^{-38}\) |

The result supports reproducibility across distinct perturbations. It does
not repair the missing screen replicate, and its empirical-null p-value must
not be described as confirmatory biological-replicate inference. In the
CRISPulator \(R=1\) benchmark this procedure is labeled BARCS-GC; ordinary
BARCS remains unchanged for replicated designs.

## Input mistakes to avoid

1. The order of `count`, `total`, and sample-table rows must match.
2. `count` cannot exceed `total`.
3. The formula is one-sided: use `~ dose + batch`, not `count ~ dose + batch`.
4. The design matrix must have full column rank.
5. Replicated FACS bins from the same donor are not independent merely because
   they appear in separate rows.
6. Multiple guides can establish perturbation consistency, but they are not
   interchangeable with independent biological screen replicates.
