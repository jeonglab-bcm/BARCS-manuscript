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
