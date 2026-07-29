# BARCS source map

The implementation follows the statistical workflow in the same order that a
user encounters it:

1. `bbreg()` validates one guide, builds the sample design matrix, and fits the
   beta-binomial mean model.
2. `bb_contrast()` tests a named linear combination of fitted coefficients.
3. `bb_screen()` applies the single-guide fit to a count matrix and returns one
   tidy row per guide.
4. `bb_calibrate_controls()` uses prespecified negative-control guides to
   estimate an empirical null scale.
5. `bb_moderate_dispersion()` stabilizes noisy guide-level dispersion
   estimates by borrowing information across the screen.
6. The `bb_gene_*()` functions provide optional, explicitly labeled
   guide-to-gene summaries.

The low-level weighted cross-products and symmetric solve are in
`src/barcs.cpp`. The R implementation remains the reference path and is used
automatically when the compiled routines are unavailable.

When reading or extending the code, keep three data contracts in view:

- `counts` is a guide-by-library matrix.
- `totals` contains the unfiltered mapped-guide total for each library and
  must not be recomputed after guide filtering.
- `data` has one row per library, in exactly the same order as the columns of
  `counts`.

The shortest executable example is
[`examples/barcs_quickstart.R`](../examples/barcs_quickstart.R).
