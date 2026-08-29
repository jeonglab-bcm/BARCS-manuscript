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
- `schmidt_tcell/BARCS_summary.md` (and the comparison tables under
  `schmidt_tcell/comparison/`)

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

The `schmidt_tcell/` directory holds a controlled head-to-head between BARCS and
the published MAGeCK pipeline on the Schmidt et al. genome-wide CRISPRa and
CRISPRi screens in primary human T cells (Science 2022,
doi:10.1126/science.abj4008), from the raw sgRNA read counts in GEO GSE174255.
Both methods are run here on the identical raw counts, so every difference is
the model rather than the input.  BARCS fits one guide-level beta-binomial
regression per library set over all twelve sorted and unsorted libraries
(`~ donor + assay * bin`), moderates the per-guide dispersions across the guides
sharing that design, and reads both cytokine contrasts off the single fit; it
never rescales the two library sets onto a common normalisation, whereas the
published pipeline normalises each library to its total, merges the two sets,
and runs paired MAGeCK RRA.  The two methods order the screens much the same way
(Spearman 0.84-0.85) and recover the paper's biology equally well - BARCS leads
inside the top 50 reference genes, MAGeCK by one gene from the top 100 outward.
They differ in calibration, measured on 660 held-out nontargeting pseudo-genes
under the same hold-out convention used elsewhere in this repository: BARCS runs
at 0.52x its nominal null rate and MAGeCK at 1.24x of the matched expectation,
so BARCS's shorter hit lists are conservatism rather than weakness.  The same
held-out check justifies the dispersion-moderation step, which adds 24% more
discoveries while lowering the null rate; it is not adopted as a general default
because on the 268-guide Liu screen it costs published hits.  Only the compact
tables under `schmidt_tcell/comparison/` are tracked; the count matrices,
guide-level tables, and MAGeCK working directories are regenerated.  See
`schmidt_tcell/README.md` for the commands.

Run the scripts in `examples/` to regenerate the full result tree.
