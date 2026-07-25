# R–browser numerical parity

`R/bbreg.R` is the reference implementation for BARCS. BARCS Web implements
the same model and inference in dependency-free JavaScript so an analysis can
run entirely in a browser.

“Same results” here means numerical and scientific equivalence, not identical
binary floating-point output. The committed parity suite evaluates the full
example screen:

- 48 guides under `~ time`, `~ time + batch`, `~ time * batch`, and `~ batch`
  (192 guide fits);
- raw and control-calibrated guide statistics;
- nine shared-effect gene results from the calibrated `~ time + batch` fit;
- identical row order, convergence status, and FDR < 0.10 decisions.

## Declared absolute tolerances

| Quantity | Maximum allowed absolute difference |
|---|---:|
| Coefficient estimate | `1e-7` |
| Standard error | `2e-7` |
| t or z statistic | `2e-5` |
| Guide p-value or FDR | `5e-8` |
| Dispersion (`rho`) | `1e-8` |
| Gene p-value or FDR | `3e-7` |

For the committed fixture, the largest observed absolute differences are:

| Quantity | Largest observed absolute difference |
|---|---:|
| Guide coefficient estimate | `1.13e-10` |
| Guide standard error | `6.36e-9` |
| Guide t statistic | `3.02e-6` |
| Guide p-value | `3.50e-8` |
| Guide FDR | `1.73e-10` |
| Dispersion (`rho`) | `3.57e-11` |
| Gene coefficient estimate | `5.29e-10` |
| Gene standard error | `5.71e-10` |
| Gene statistic | `1.39e-6` |
| Gene p-value | `1.34e-8` |
| Gene FDR | `6.00e-8` |

The declared bounds cover these measured errors without hiding a
scientifically meaningful discrepancy. Values extremely close to a reporting
threshold should still be reviewed at full precision.

Small discrepancies are expected because R and JavaScript can differ in IRLS
stopping paths, linear-algebra operation order, one-dimensional root finding,
and implementations of Student-t and normal probabilities. All calculations
use IEEE-754 double precision, but equal equations do not imply bit-for-bit
equality across independent numerical implementations.

## Reproduce the reference and test

From the repository root:

```sh
Rscript web/scripts/generate-r-parity-reference.R guides \
  > web/tests/fixtures/r-guide-reference.csv
Rscript web/scripts/generate-r-parity-reference.R genes \
  > web/tests/fixtures/r-gene-reference.csv
cd web
npm test
```

The fixtures are committed so the web suite remains runnable on systems
without R. Regenerating them requires the same `R/bbreg.R` source used by the
package tests.

## Manuscript-scale parity

The GSE70038 preset is a second, independent parity layer. It uses the
manuscript's 64,747-guide count matrix, 16-library Table 5 design, no
minimum-count filter, and its explicitly labeled median-effect/directional
Stouffer comparison summary for all four terminal coefficients.

The test requires:

- 64,747 guide results and 18,077 gene results for each coefficient;
- exactly 3,608, 5,702, 7,080, and 5,737 gene discoveries at FDR < 0.10;
- the same coefficient-specific top 20 genes;
- matching coefficient-specific effect, p-value, and FDR checksums;
- matching effects, p-values, FDRs, and decisions for FBXO42, HDAC2, HEATR1,
  PKMYT1, TFAP2C, and WEE1 under every coefficient.

Across all four sets of 64,747 guide rows, the largest observed absolute
differences are `1.40e-11` for the coefficient, `5.67e-14` for the standard
error, `8.42e-12` for the p-value and FDR, and `3.05e-18` for dispersion.
Across all four sets of 18,077 genes, the largest observed differences are
`2.04e-12` for the median effect, `5.14e-9` for the p-value, and `5.26e-9` for
FDR. All FDR < 0.10 decisions are identical.

Achieving this agreement requires matching R's Brent-style root finder,
binomial `glm.fit` initialization, stable deviance arithmetic, and probability
tails. These are numerical implementation details, not changes to the BARCS
estimating equations.
