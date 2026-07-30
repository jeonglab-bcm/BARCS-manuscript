# Generated benchmark results

Benchmark scripts write their complete outputs under this directory. Large
count matrices, per-guide tables, and third-party intermediate files are
excluded from Git history.

The repository tracks compact summary tables used by the manuscript:

- `simulation_summary.csv`
- `simulation_diagnostics.csv`
- `gse70038/method_concordance.csv`
- `gse70038/published_gene_results.csv`
- `chronos_tzelepis/benchmark_metrics.csv`
- `chronos_tzelepis/cnv_bias_diagnostic.csv`
- `chronos_tzelepis/effect_rank_correlations.csv`
- `sanson_benchmark/benchmark_metrics.csv`
- `sanson_benchmark/cnv_bias_diagnostic.csv`
- `sanson_benchmark/paired_bootstrap_comparison.csv`
- `waterbear_facs/benchmark_metrics.csv`
- `waterbear_facs/effect_rank_correlations.csv`
- `waterbear_facs/non_targeting_calibration.csv`
- `waterbear_facs/validation_panel_metrics.csv`

The compact simulation diagnostics report the guide-level null rejection
rate and the fraction of fitted dispersions at the lower boundary.  The
manuscript's Waterbear null-calibration result is the deterministic five-fold
held-out estimate; a seeded fold permutation is retained as a structure check.
Fitting a scale and evaluating it on the same controls is not treated as
validation.

The Liang benchmark applies one common signed-normal tail-scaling rule to all
five methods after their native gene summaries.  Non-targeting guides are
assigned to the same deterministic pseudo-genes for every method, with
pseudo-gene sizes sampled from the observed target-gene guide-count
distribution.  The scale therefore includes the native guide-to-gene
aggregation step.  The accompanying time-point ablation compares BARCS fits
using days 0, 7, and 14 with fits using days 0 and 14 only.

The null calibration grid varies the number of independent libraries, guide
abundance, overdispersion through the continuous-dose failure regime, and
three versus five guides per gene.  A focused arm adds within-gene guide
correlation, and a deterministic null-gene split measures post-aggregation
held-out calibration.  Seed-level and summarized guide- and gene-error rates
are tracked under `data/derived/`, together with the observed lower-boundary
fractions for HT-29, IL2RA, and Liang.

Run the scripts in `examples/` to regenerate the full result tree.
