# CB² to CB²-Reg: a computational walkthrough

This version proves the relationship by running small examples. You can execute
everything at once from the repository root:

```sh
Rscript examples/cb2_generalization_walkthrough.R
```

The complete runnable file is
[`examples/cb2_generalization_walkthrough.R`](../examples/cb2_generalization_walkthrough.R).

## What we are checking

Original CB² calculates

```text
effect = weighted proportion in B - weighted proportion in A

                     effect
t = -----------------------------------------
    sqrt(variance of A + variance of B)
```

We want to check three claims:

1. A weighted regression with a 0/1 group column gives exactly the same effect,
   standard error, and t-statistic.
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

The direct CB² calculation is:

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

This is the entire exact proof in one small computation.

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

The exact result above embeds the original CB² summaries in weighted
regression. The default `bbreg()` function makes additional modeling choices:

- logit effect instead of raw proportion difference;
- one guide-level dispersion across the design instead of separate group fits;
- residual degrees of freedom instead of the original group-variance formula.

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

```r
p0 <- 0.002
delta <- c(8e-4, 4e-4, 2e-4, 1e-4, 5e-5, 2.5e-5)
```

The walkthrough compares the raw-proportion t-statistic with its delta-method
logit version:

| B − A | Raw t | Logit t | Logit/raw | Absolute gap |
|---:|---:|---:|---:|---:|
| 0.000800 | 4.697762 | 4.639446 | 0.987587 | 0.058316 |
| 0.000400 | 2.348881 | 2.362284 | 1.005706 | 0.013403 |
| 0.000200 | 1.174440 | 1.181210 | 1.005764 | 0.006770 |
| 0.000100 | 0.587220 | 0.589332 | 1.003597 | 0.002112 |
| 0.000050 | 0.293610 | 0.294190 | 1.001975 | 0.000580 |
| 0.000025 | 0.146805 | 0.146957 | 1.001032 | 0.000151 |

As the group difference approaches zero:

```text
logit t / raw t  approaches 1
absolute gap     approaches 0
```

This is the computational version of the local-equivalence proof.

## Example 6: the part that is genuinely new

In the original two-group problem, the second design-matrix column is:

```text
0, 0, 0, 0, 1, 1, 1, 1
```

That column only asks whether B differs from A.

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

That is the practical meaning of generalization.

## What each example establishes

| Example | What it shows |
|---|---|
| 1 | The direct two-group calculation and weighted regression are identical |
| 2 | The identity survives 1,000 random numerical stress tests |
| 3 | The identity holds after original `fit_ab()` estimates from read counts |
| 4 | Default logit `bbreg()` is close but not exactly equal in finite samples |
| 5 | Raw and logit statistics converge as the group difference approaches zero |
| 6 | The 0/1 group column can be replaced by continuous dose |

Computational checks cannot replace the algebraic proof for every possible
dataset, but they make each step observable. The formal proof is in
[`main.tex`](../main.tex), and the simpler conceptual explanation is in
[`cb2-generalization-high-school-proof.md`](cb2-generalization-high-school-proof.md).
