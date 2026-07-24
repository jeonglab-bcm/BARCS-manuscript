# CB² to CB²-Reg: a computational walkthrough

This version demonstrates the relationship with small examples. Computational
checks verify the code paths; the manuscript supplies the algebra and
asymptotic assumptions. You can execute everything at once from the repository
root:

```sh
Rscript examples/cb2_generalization_walkthrough.R
```

The complete runnable file is
[`examples/cb2_generalization_walkthrough.R`](../examples/cb2_generalization_walkthrough.R).

## What we are checking

Original CB² calculates

$$
\widehat{\Delta}_{\mathrm{CB2}}
  = \widehat p_B-\widehat p_A,
\qquad
T_{\mathrm{CB2}}
  = \frac{\widehat p_B-\widehat p_A}
  {\sqrt{V_A+V_B}}.
$$

Here, $\widehat p_A$ and $\widehat p_B$ are the weighted guide proportions, and
$V_A$ and $V_B$ are their estimated variances.

We want to check three claims:

1. After original CB² has computed its group summaries, a saturated weighted
   regression with a 0/1 group column reproduces the same effect, standard
   error, and t-statistic.
2. The default logit regression is close near the null hypothesis, although it
   need not be exactly equal in a small dataset.
3. Once the group column has been written as a design-matrix column, we can
   replace 0/1 with dose, time, or another quantitative variable.

## Example 1: numbers that can be checked by hand

Suppose original CB² has already estimated:

```r
p <- c(A = 0.0014, B = 0.0023)
v <- c(A = 4e-8, B = 5e-8)
```

In mathematical notation, the inputs are

$$
\widehat p_A=0.0014,\quad
\widehat p_B=0.0023,\quad
V_A=4\times10^{-8},\quad
V_B=5\times10^{-8}.
$$

The direct CB² calculation is:

$$
\begin{aligned}
\widehat\Delta
  &= 0.0023-0.0014=0.0009,\\
\mathrm{SE}(\widehat\Delta)
  &= \sqrt{4\times10^{-8}+5\times10^{-8}}
   =0.0003,\\
T_{\mathrm{CB2}}
  &=\frac{0.0009}{0.0003}=3.
\end{aligned}
$$

```r
effect <- p["B"] - p["A"]
se <- sqrt(v["A"] + v["B"])
t <- effect / se

c(effect = effect, standard_error = se, t_value = t)
```

Output:

```text
        effect standard_error        t_value
        0.0009         0.0003          3.0000
```

Now create the smallest possible two-group design matrix:

```r
x <- cbind(
  intercept = 1,
  group_B = c(0, 1)
)
x
```

```text
     intercept group_B
[1,]         1       0
[2,]         1       1
```

The first row represents group A, and the second represents group B. Give the
two rows precision weights `1 / VA` and `1 / VB`:

$$
X=
\begin{pmatrix}
1&0\\
1&1
\end{pmatrix},
\qquad
\mathbf y=
\begin{pmatrix}
\widehat p_A\\
\widehat p_B
\end{pmatrix},
\qquad
\Omega=
\begin{pmatrix}
1/V_A&0\\
0&1/V_B
\end{pmatrix}.
$$

The weighted-regression calculation is

$$
\widehat{\boldsymbol\beta}
  =(X^\top\Omega X)^{-1}X^\top\Omega\mathbf y,
\qquad
\widehat{\mathrm{Cov}}(\widehat{\boldsymbol\beta})
  =(X^\top\Omega X)^{-1}.
$$

```r
precision <- 1 / v
information <- crossprod(x, precision * x)
score_target <- crossprod(x, precision * p)

beta <- solve(information, score_target)
covariance <- solve(information)

beta
sqrt(covariance[2, 2])
beta[2] / sqrt(covariance[2, 2])
```

Output:

```text
intercept   group_B
   0.0014    0.0009

standard error of group_B = 0.0003
t-statistic               = 3
```

The regression intercept is the group-A proportion. The regression slope is
`B - A`. Its standard error and t-statistic are exactly the direct CB² values.

In symbols, the computation has produced

$$
\widehat{\boldsymbol\beta} =
\begin{pmatrix}
\widehat p_A\\
\widehat p_B-\widehat p_A
\end{pmatrix},
\qquad
\widehat{\mathrm{Cov}}(\widehat{\boldsymbol\beta}) =
\begin{pmatrix}
V_A&-V_A\\
-V_A&V_A+V_B
\end{pmatrix}.
$$

Therefore,

$$
\frac{\widehat\beta_1}
{\sqrt{\widehat{\mathrm{Var}}(\widehat\beta_1)}} =
\frac{\widehat p_B-\widehat p_A}{\sqrt{V_A+V_B}}
=T_{\mathrm{CB2}}.
$$

This is the legacy representation lemma in one small computation. It is a
summary-scale compatibility identity, not proof that default `bbreg()` is a
strictly nested version of the original estimator.

## Example 2: try to break the identity

One successful example might be a coincidence, so the walkthrough generates
1,000 random pairs of proportions and variances:

```r
set.seed(20260723)

random_error <- replicate(1000, {
  random_p <- runif(2, 1e-5, 0.05)
  variance_scale <- 10^runif(1, -10, -6)
  random_v <- variance_scale * runif(2, 0.5, 2)
  random_precision <- 1 / random_v

  information <- crossprod(x, random_precision * x)
  beta <- solve(
    information,
    crossprod(x, random_precision * random_p)
  )
  covariance <- solve(information)

  direct_effect <- random_p[2] - random_p[1]
  direct_se <- sqrt(sum(random_v))

  max(abs(c(
    beta[1] - random_p[1],
    beta[2] - direct_effect,
    covariance[2, 2] - sum(random_v),
    beta[2] / sqrt(covariance[2, 2]) -
      direct_effect / direct_se
  )))
})

max(random_error)
```

Every random problem uses a different abundance and uncertainty. The reported
maximum discrepancy is at floating-point rounding level. The formal algebraic
proof covers every positive variance; this randomized check shows that the
code actually performs that algebra.

For random problem $r$, the program records the largest of four discrepancies:

$$
\epsilon_r=\max\left\{
\left|\widehat\beta_{0r}-\widehat p_{Ar}\right|,
\left|\widehat\beta_{1r}-(\widehat p_{Br}-\widehat p_{Ar})\right|,
\left|\widehat{\mathrm{Var}}(\widehat\beta_{1r})-(V_{Ar}+V_{Br})\right|,
\left|T_{\mathrm{GLS},r}-T_{\mathrm{CB2},r}\right|
\right\}.
$$

The final check is $\max_{1\le r\le1000}\epsilon_r$, which is approximately
$1.02\times10^{-12}$ on the recorded run.

## Example 3: start from read counts and run original CB²

Now use four samples in each group:

```r
count_a <- matrix(c(74, 112, 91, 139), nrow = 1)
total_a <- matrix(c(52000, 81000, 69000, 97000), nrow = 1)

count_b <- matrix(c(128, 177, 164, 211), nrow = 1)
total_b <- matrix(c(55000, 76000, 72000, 93000), nrow = 1)

fit_a <- CB2::fit_ab(count_a, total_a)
fit_b <- CB2::fit_ab(count_b, total_b)

phat <- c(A = fit_a$phat, B = fit_b$phat)
vhat <- c(A = fit_a$vhat, B = fit_b$vhat)
```

For sample $i$, original CB² first forms

$$
Y_i=\frac{K_i}{N_i}.
$$

Within group $g$, `fit_ab()` estimates normalized weights $w_{gi}$ and then
computes

$$
\widehat p_g=\sum_{i\in g}w_{gi}Y_i,
\qquad
\sum_{i\in g}w_{gi}=1,
$$

together with the group variance $V_g$. Thus the inputs to the regression in
this example are estimated from the read counts rather than typed by hand.

Original CB² converts the shared t-statistic to a p-value with

$$
\nu_{\mathrm{CB2}} =
\frac{(V_A+V_B)^2}
{V_A^2/(n_A-1)+V_B^2/(n_B-1)}
$$

degrees of freedom and

$$
p_{\text{two-sided}}
=2\,\Pr\!\left(t_{\nu_{\mathrm{CB2}}}
\geq\left|T_{\mathrm{CB2}}\right|\right).
$$

Using these same degrees of freedom for the regression slope therefore makes
the p-values identical as well.

The walkthrough feeds those original CB² estimates into the same two-row
weighted regression. The result is:

| Quantity | Original CB² | Weighted regression |
|---|---:|---:|
| Effect | 0.000905992949471 | 0.000905992949471 |
| Standard error | 0.000111310501148 | 0.000111310501148 |
| t-statistic | 8.139330432681405 | 8.139330432681406 |
| Degrees of freedom | 5.647255017373669 | 5.647255017373669 |
| Two-sided p-value | 0.000251292755657 | 0.000251292755657 |

The maximum difference is approximately `1.78e-15`, which is computer rounding
error.

This example is stronger than starting with invented group proportions:
`fit_ab()` first estimates the weighted proportions and beta-binomial
variances from actual read counts.

## Example 4: deliberately compare the two finite-sample estimators

The representation above embeds the completed original CB² summaries in
weighted regression. The default `bbreg()` function makes additional modeling
choices:

- logit effect instead of raw proportion difference;
- one guide-level dispersion across the design instead of separate group fits;
- residual degrees of freedom instead of the original group-variance formula.

The two effect parameters are

$$
\Delta_{\mathrm{CB2}}=p_B-p_A
$$

and

$$
\beta_{\mathrm{group}}
  =\mathrm{logit}(\mu_B)-\mathrm{logit}(\mu_A)
  =\log\left(\frac{\mu_B/(1-\mu_B)}
  {\mu_A/(1-\mu_A)}\right).
$$

They answer the same directional question but use different units.

Run:

```r
sample_data <- data.frame(
  group = factor(rep(c("A", "B"), each = 4))
)

logit_fit <- bbreg(
  count = c(count_a, count_b),
  total = c(total_a, total_b),
  formula = ~ group,
  data = sample_data
)

logit_fit$coefficient_table["groupB", ]
```

For these data:

| Method | Effect | Effect scale | t | df | p |
|---|---:|---|---:|---:|---:|
| Original CB² | 0.000905993 | Raw proportion | 8.13933 | 5.64726 | 0.000251293 |
| Default CB²-Reg | 0.502399 | Log odds | 8.06432 | 6 | 0.000194616 |

The effect numbers should not be compared directly because they have different
units. The test statistics are close but not identical. This is expected and
is why the manuscript does **not** claim that `bbreg(~ group)` reproduces
`measure_sgrna_stats()` exactly in every finite dataset.

## Example 5: move toward the null hypothesis

Set a common guide proportion of `0.002`, then place groups A and B equally
around it. Repeatedly reduce the distance between them:

$$
p_A=p_0-\frac{\delta}{2},
\qquad
p_B=p_0+\frac{\delta}{2},
\qquad p_0=0.002.
$$

```r
p0 <- 0.002
delta <- c(8e-4, 4e-4, 2e-4, 1e-4, 5e-5, 2.5e-5)
```

The walkthrough compares the raw-proportion t-statistic with its delta-method
logit version:

$$
T_{\mathrm{raw}}
=\frac{p_B-p_A}{\sqrt{V_A+V_B}},
$$

$$
T_{\logit} =
\frac{\mathrm{logit}(p_B)-\mathrm{logit}(p_A)}
{\sqrt{\{\mathrm{logit}'(p_A)\}^2V_A+
\{\mathrm{logit}'(p_B)\}^2V_B}},
\qquad
\mathrm{logit}'(p)=\frac{1}{p(1-p)}.
$$

| B − A | Raw t | Logit t | Logit/raw | Absolute gap |
|---:|---:|---:|---:|---:|
| 0.000800 | 4.697762 | 4.639446 | 0.987587 | 0.058316 |
| 0.000400 | 2.348881 | 2.362284 | 1.005706 | 0.013403 |
| 0.000200 | 1.174440 | 1.181210 | 1.005764 | 0.006770 |
| 0.000100 | 0.587220 | 0.589332 | 1.003597 | 0.002112 |
| 0.000050 | 0.293610 | 0.294190 | 1.001975 | 0.000580 |
| 0.000025 | 0.146805 | 0.146957 | 1.001032 | 0.000151 |

As the group difference approaches zero and the number of independent
biological libraries grows:

$$
\frac{T_{\logit}}{T_{\mathrm{raw}}}\longrightarrow1,
\qquad
\left|T_{\logit}-T_{\mathrm{raw}}\right|\longrightarrow0.
$$

The reason is the first-order Taylor approximation

$$
\mathrm{logit}(p_B)-\mathrm{logit}(p_A)
=\mathrm{logit}'(p_0)(p_B-p_A)+o(\delta).
$$

The standard error is multiplied by the same
$\mathrm{logit}'(p_0)$, so that factor
cancels from the t-statistic.

This is the computational illustration of the local-equivalence result.
Increasing sequencing depth while holding the number of independent libraries
fixed is not the asymptotic regime used by that result.

## Example 6: the part that is genuinely new

In the original two-group problem, the second design-matrix column is:

```text
0, 0, 0, 0, 1, 1, 1, 1
```

That column only asks whether B differs from A.

The corresponding mean model is

$$
\mathrm{logit}(\mu_i)=\beta_0+\beta_1G_i,
\qquad G_i\in\{0,1\}.
$$

For a dose experiment, replace it with:

```text
0, 0, 1, 1, 2, 2, 3, 3, 4, 4
```

The design matrix becomes:

```r
dose <- rep(0:4, each = 2)
model.matrix(~ dose)
```

```text
      (Intercept) dose
 [1,]           1    0
 [2,]           1    0
 [3,]           1    1
 [4,]           1    1
 [5,]           1    2
 [6,]           1    2
 [7,]           1    3
 [8,]           1    3
 [9,]           1    4
[10,]           1    4
```

The new mean model is

$$
\mathrm{logit}(\mu_i)=\beta_0+\beta_{\mathrm{dose}}d_i,
\qquad d_i\in\{0,1,2,3,4\}.
$$

More generally, CB²-Reg uses

$$
\mathrm{logit}(\mu_i)=\mathbf x_i^\top\boldsymbol\beta,
$$

where $\mathbf x_i$ can also contain batch, donor, CNV, interaction, or spline
columns.

The example creates counts with a true logit slope of `0.35` and fits:

```r
dose_fit <- bbreg(
  count = dose_count,
  total = dose_total,
  formula = ~ dose,
  data = data.frame(dose = dose)
)
```

Result:

```text
             estimate  standard error  t value  df       p
(Intercept) -7.000958       0.052623  -133.04   8 1.14e-14
dose         0.350251       0.017729    19.76   8 4.49e-08
```

The fitted dose effect is `0.350251`, very close to the value `0.35` used to
create the data.

Nothing fundamental changed in the data table:

- the first column still represents the baseline;
- the second column now represents dose instead of a 0/1 group label;
- more columns can represent batch, donor, CNV, or interactions.

That is the practical meaning of the design-matrix extension.

## What each example establishes

| Example | What it shows |
|---|---|
| 1 | The completed legacy summaries have an identical saturated GLS representation |
| 2 | The identity survives 1,000 random numerical stress tests |
| 3 | The identity holds after original `fit_ab()` estimates from read counts |
| 4 | Default logit `bbreg()` is close but not exactly equal in finite samples |
| 5 | Raw and logit statistics converge as the group difference approaches zero |
| 6 | The 0/1 group column can be replaced by continuous dose |

Computational checks cannot replace the algebraic proof for every possible
dataset, but they make each step observable. The formal proof is in
[`main.tex`](../main.tex), and the simpler conceptual explanation is in
[`cb2-generalization-high-school-proof.md`](cb2-generalization-high-school-proof.md).

The default covariance in these examples treats the estimated guide
dispersion as a fixed plug-in value. Its Student t reference does not formally
propagate dispersion-estimation uncertainty; this is a separate small-sample
limitation from the legacy-representation question.
