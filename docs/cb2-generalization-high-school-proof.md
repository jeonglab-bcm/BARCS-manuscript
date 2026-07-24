# How original CB² connects to CB²-Reg

If you prefer to learn by running code, use the
[`computational walkthrough`](cb2-generalization-computational-walkthrough.md).
It starts with a three-line calculation, runs 1,000 random stress tests, calls
the original `fit_ab()`, and finishes with a continuous-dose example.

## The short answer

Original CB² compares two groups. CB²-Reg can analyze two groups, but it can
also analyze dose, time, batch, donor, and other variables together.

For the simple two-group problem, the completed old CB² group summaries can be
written as a saturated two-column regression calculation. The effect, standard
error, t-statistic, and p-value are then exactly the same. This is a useful
compatibility representation, but it does not prove that the default
`bbreg()` estimator is literally the old estimator with extra columns.

The default CB²-Reg model uses a logit transformation, so its answer is not
always numerically identical to old CB² in a small dataset. However, when the
two groups are close—the situation considered by the null hypothesis—the two
test statistics approach each other.

That gives us two different conclusions:

1. **Legacy representation:** after old CB² has computed its weighted group
   proportions and variances, a two-cell weighted regression reproduces its
   statistic exactly.
2. **Local-equivalence result:** the default logit version of CB²-Reg
   approaches the old CB² test for small group differences as the number of
   independent biological libraries grows.

## Step 1: What original CB² does

Imagine that we have a control group called A and a treatment group called B.
For one guide, each sample gives us a proportion:

```text
guide reads / all mapped reads
```

Some samples are sequenced more deeply or are less noisy than others. CB²
therefore does not treat every sample equally. It calculates a weighted average
for each group:

```text
pA = weighted average guide proportion in group A
pB = weighted average guide proportion in group B
```

The weights within each group add to 1. CB² also estimates the uncertainty of
the two weighted averages:

```text
VA = variance of pA
VB = variance of pB
```

The estimated effect is

```text
effect = pB - pA
```

Because groups A and B are independent, the variance of their difference is

```text
variance of effect = VA + VB
```

Therefore, the standard error is

```text
standard error = square root of (VA + VB)
```

and the original CB² statistic is

```text
                 pB - pA
t_CB2 = ---------------------------
          square root of (VA + VB)
```

## Step 2: Turn the group labels into two regression columns

A design matrix is just a table of numbers that describes the samples. For two
groups, it can be this simple:

| Sample group | Intercept column | Group-B column |
|---|---:|---:|
| A | 1 | 0 |
| A | 1 | 0 |
| B | 1 | 1 |
| B | 1 | 1 |

The intercept column is always 1. The Group-B column is 0 for group A and 1
for group B.

The regression equation is

```text
predicted proportion = intercept + slope × Group-B
```

For group A, `Group-B = 0`, so

```text
predicted A proportion = intercept
```

For group B, `Group-B = 1`, so

```text
predicted B proportion = intercept + slope
```

To match the two CB² weighted averages, the regression must therefore choose

```text
intercept = pA
slope     = pB - pA
```

This is already the main idea of the proof: the old CB² group difference is
the slope of a two-group regression.

## Step 3: Give the regression the CB² uncertainties

We give group A total statistical weight `1 / VA` and group B total
statistical weight `1 / VB`. Within each group, we divide that total weight
using the same sample weights that CB² calculated.

This makes the regression remember both parts of the old calculation:

- the weighted average in each group;
- the uncertainty of that weighted average.

Solving this weighted regression gives

```text
regression intercept = pA
regression slope     = pB - pA
variance of slope    = VA + VB
```

Therefore,

```text
                         regression slope
t_regression = -------------------------------------
                square root of variance of the slope

                         pB - pA
             = ---------------------------
                square root of (VA + VB)

             = t_CB2
```

If we also use the same CB² degrees-of-freedom formula, the t-distribution is
the same. Therefore, the p-value is exactly the same.

This shows that original CB² has an exact weighted-regression representation
on its group-summary scale. Because any two independent group summaries can be
placed in a saturated two-cell regression, this calculation alone does not
prove that default `bbreg()` strictly contains the old estimator.

## Step 4: Why the default logit model is slightly different

A proportion must stay between 0 and 1. A straight-line regression can
sometimes predict impossible values such as `-0.01` or `1.04`. CB²-Reg avoids
this problem by using a curved transformation called the logit.

The important high-school-calculus idea is that a smooth curve looks almost
like a straight line when we zoom in far enough.

Near a common proportion `p`, the logit transformation changes a small
difference approximately like this:

```text
logit(pB) - logit(pA)

                         1
approximately = ------------------- × (pB - pA)
                    p × (1 - p)
```

The logit standard error is multiplied by approximately the same number:

```text
                         1
logit standard error ≈ ----------- × square root of (VA + VB)
                    p × (1 - p)
```

When we form the t-statistic, that number appears on both the top and bottom:

```text
       [1 / {p(1-p)}] × (pB - pA)
t ≈ ---------------------------------
       [1 / {p(1-p)}] × sqrt(VA + VB)
```

The common multiplier cancels:

```text
             pB - pA
t ≈ ---------------------------
      square root of (VA + VB)
```

That is the original CB² statistic.

So the logit version is not exactly equal for every small dataset, but the two
statistics become equivalent as the group difference becomes small and the
number of independent biological libraries grows. Sequencing the same fixed
set of libraries more deeply is not enough for this asymptotic argument.

## Step 5: Why the sample weights also match

Let `N` be a sample's sequencing depth and let `kappa` describe how much
between-sample variation exists.

The old CB² relative weight has the form

```text
kappa × N
-----------
kappa + N
```

The CB²-Reg relative weight has the form

```text
    N
-----------
kappa + N
```

The only difference is the common `kappa` multiplier. Multiplying every weight
in one group by the same number does not change a weighted average because the
weights are normalized to add to 1.

For example:

```text
weights:           2, 4, 6
normalized:      1/6, 2/6, 3/6

multiply by 10:  20, 40, 60
normalized:      1/6, 2/6, 3/6
```

Therefore, when the old and new models estimate the same underlying
dispersion, they use the same relative depth-aware weighting.

## A numerical example

The executable example in
[`examples/cb2_generalization_proof.R`](../examples/cb2_generalization_proof.R)
uses four samples in each group.

It obtains:

| Quantity | Original CB² | Weighted regression |
|---|---:|---:|
| Effect | 0.000905992949471 | 0.000905992949471 |
| Standard error | 0.000111310501148 | 0.000111310501148 |
| t-statistic | 8.139330432681405 | 8.139330432681406 |
| Degrees of freedom | 5.647255017373669 | 5.647255017373669 |
| Two-sided p-value | 0.000251292755657 | 0.000251292755657 |

The largest difference between the two calculations is

```text
0.00000000000000178
```

This is rounding error from the computer, not a statistical difference.

The same example moves two proportions closer together. As their difference
shrinks from `0.0008` to `0.000025`, the gap between the raw-proportion and
logit t-statistics shrinks from `0.0583` to `0.000151`.

## What has and has not been proved

### What is established

- The completed original CB² group summaries have an exact saturated
  weighted-regression representation.
- Its standard error and p-value are exactly recovered when we use its original
  group variances and degrees of freedom.
- Under a common beta-binomial dispersion, old CB² and CB²-Reg use the same
  relative depth-aware weights.
- The default logit statistic approaches the original statistic for local,
  small group differences.
- A design matrix lets the same regression idea represent dose, time, batch,
  donor, interactions, and adjusted contrasts.

### What is not claimed

- The representation lemma does not prove strict nesting of the implemented
  `bbreg()` estimator.
- `bbreg(~ group)` is not guaranteed to return exactly the same finite-sample
  number as `measure_sgrna_stats()`.
- Original CB² estimates the two groups separately, while CB²-Reg estimates one
  guide-level dispersion across the entire design.
- Original CB² uses a group-variance degrees-of-freedom formula, while CB²-Reg
  normally uses residual sample degrees of freedom.
- The model-based CB²-Reg covariance does not propagate uncertainty in the
  per-guide dispersion estimate. Student t degrees of freedom are not a formal
  correction for that uncertainty.
- Greater read depth cannot replace additional independent biological
  libraries in the local-equivalence argument.
- Neither proof says CB²-Reg must outperform MAGeCK, Waterbear, MAUDE, Chronos,
  or another specialized method on every dataset.

## The one-sentence conclusion

Original CB² has a compatible two-cell regression representation, while
CB²-Reg extends the beta-binomial weighting principle to coefficient-based
questions; under explicit common-dispersion and replication assumptions, its
two-group logit test agrees with CB² to first order near the null hypothesis.

For the formal matrix derivation and assumptions, see
[`main.tex`](../main.tex) or the
[`rendered manuscript`](../output/pdf/beta-binomial-regression-continuous-phenotypes.pdf).
