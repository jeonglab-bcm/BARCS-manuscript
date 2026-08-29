# BARCS versus the published MAGeCK analysis of the Schmidt et al. T cell CRISPRa/CRISPRi screens

**Prepared by:** Hyun-Hwan Jeong (hyun-hwan.jeong@bcm.edu)
**Date:** 29 August 2026
**Source:** Schmidt R, Steinhart Z, Layeghi M, *et al.* "CRISPR activation and interference screens decode stimulation responses in primary human T cells." *Science* 376, eabj4008 (2022). doi:10.1126/science.abj4008 — raw sgRNA read counts from GEO [GSE174255](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174255).

---

## One-paragraph summary

This is a controlled head-to-head, not a re-analysis: BARCS and the paper's own
MAGeCK pipeline are both run here, on the same GEO raw counts, so every
difference is the model rather than the input. BARCS fits one guide-level
beta-binomial regression per library set across all twelve sorted and unsorted
libraries, with donor as an explicit term, moderates the per-guide dispersions
across the ~113,000 guides sharing that design, and reads both cytokine
contrasts off the single fit; the published pipeline normalises each library to
its total, merges the two library sets, and runs a paired MAGeCK RRA test per
screen. The fits are clean — 99.9% guide convergence, and 0 of 165 nontargeting
pseudo-genes called in every one of the four screens — and BARCS reproduces the
paper's biology: IFNG and IL2 top their own screens, VAV1/CD28/LCP2 lead the
positive regulators and MAP4K1/SLA2/GRAP the negative ones, TBX21 is up and
GATA3 down in IFN-γ only, and the MALT1–BCL10–TRAF6–MAP3K7–IκB circuit appears
in CRISPRi IFN-γ and nowhere else, exactly as published. **On recovering that
biology the two methods are level**: BARCS leads inside the top 50 (36 versus 34
of 43 reference assertions), MAGeCK by one from the top 100 outward. Where they
differ is calibration. On 660 held-out nontargeting pseudo-genes BARCS's null
p-values run at 0.52× their nominal rate and MAGeCK's at 1.24× of the matched
expectation, so BARCS's shorter hit lists are conservatism rather than weakness,
and being at half nominal means power is still being left on the table.

## What was done

| Item | Description |
|---|---|
| Data | GEO GSE174255, `GSE174255_sgRNA-Read-Counts.xlsx` — raw integer per-library sgRNA counts, the input the authors ran MAGeCK on |
| CRISPRa | Calabrese Set A (56,762 guides, 54.7M reads) + Set B (56,476, 47.5M); 18,805 genes analysed, 3 + 3 guides each |
| CRISPRi | Dolcetto Set A (57,050 guides, 61.1M reads) + Set B (57,011, 52.8M); 18,865 genes analysed |
| Controls | 496 nontargeting sgRNAs per library set, disjoint between sets — 992 per screen |
| Libraries | 12 per set: 2 donors × {IL-2, IFN-γ} × {low, high, unsorted}. The plasmid library is pre-transduction and is excluded |
| Readout | CRISPRa/CRISPRi T cells sorted on IL-2 (CD4⁺) or IFN-γ (CD8⁺), high versus low bins |
| BARCS | `bb_screen(~ donor + assay * bin)` per library set → `bb_moderate_dispersion()` → `bb_calibrate_controls()` → `bb_gene_consistency()` |
| MAGeCK | The Methods recipe: normalise to library total, merge both sets, `mageck test --norm-method none --paired --control-sgrna` |

A positive effect means the perturbation enriched the guide in the
cytokine-high bin, matching the paper's log2(high/low) sign convention.

**The two designs, and why they differ.** The published pipeline has to put two
separately sequenced library pools on a common scale before it can compare
anything: it divides every library by its own total, merges the Set A and Set B
tables, then tells MAGeCK not to normalise again. BARCS never needs that step. A
guide's high-versus-low effect is a contrast internal to the libraries that guide
was sequenced in, taken against those libraries' own read totals, so Set A and
Set B guides are directly comparable without either set being rescaled onto the
other; the sets meet only at the gene level, where they contribute 3 + 3 = 6
guides. Within a set, one model per guide covers all twelve libraries, which
keeps donor in the design rather than using it as a pairing device, lets the
unsorted bins anchor the dispersion, and yields both cytokine contrasts from a
single fit. `bb_screen()` reports one model-matrix coefficient, so the fit is run
once per cytokine with `assay` releveled; the two runs are the same model in two
parameterisations.

**Dispersion moderation is part of the pipeline, and is not a free default.** A
per-guide dispersion estimated on 5 residual degrees of freedom is badly
determined. `bb_moderate_dispersion()` borrows across the guides sharing the
design, raising the effective df from 5 to ~140 (CRISPRa) and ~107 (CRISPRi).
The held-out check in Result 4 shows it raises discoveries from 1,677 to 2,082
while *lowering* the held-out null rate from 0.036 to 0.026 — power and
calibration together, not a trade. It is applied here because this screen has the
guide count to support it, and deliberately not adopted as a general default:
on the 268-guide Liu screen in `results/liu_tcell/` the same step costs four
published hits (TNFAIP3, NFKBIA, GNAS, RASA2) and *raises* the null scale from
1.72 to 1.95, because the prior is estimated from too few guides.

**The MAGeCK arm is a faithful reproduction.** Run on the GEO counts with the
published options it returns 480, 563, 238 and 211 hits under the paper's own
rule (median |log2FC| > 0.5 and FDR < 0.05) against the 444, 471, 226 and 203
reported in the text — within 8–19%, the residual attributable to MAGeCK 0.5.9.5
here versus 0.5.9.2 there and to the one detail the Methods leave implicit (see
Limits #2). Counts quoted below are 477, 562, 238 and 211, the same runs
restricted to genes both methods score.

## Result 1 — the fits are trustworthy

| Check | CRISPRa (IL-2 / IFN-γ) | CRISPRi (IL-2 / IFN-γ) |
|---|---|---|
| Guides fitted | 113,238 | 114,061 |
| IRLS convergence | 99.90% | 99.92% |
| Residual df per guide, before → after moderation | 5 → 140 | 5 → 107 |
| Median extra-binomial ρ | 3.7 × 10⁻⁷ | 1.6 × 10⁻⁷ |
| Guide-level control calibration | × 1.43 / × 1.41 | × 1.05 / × 1.05 |
| Gene-level null scale | 1.60 / 1.68 | 1.15 / 1.19 |
| Nontargeting pseudo-genes at FDR < 0.05 | 0 / 165 · 0 / 165 | 0 / 165 · 0 / 165 |

The controls are grouped into 165 pseudo-genes shaped like real genes — three
Set A sgRNAs and three Set B sgRNAs, six guides — rather than left as the single
992-guide `NO-TARGET` gene the published analysis uses. One 992-guide gene cannot
report how often a *gene-sized* group of null guides is called, and it leaves
`bb_gene_consistency()` with one control gene, too few to estimate a
control-based null at all. Not one of the 660 control pseudo-genes across the
four screens is called.

## Result 2 — the published biology reproduces

Top of each screen, unlabelled and unfiltered:

- **CRISPRa IFN-γ**: IFNG (rank 1, the screen's own positive control), VAV1, CD28,
  TNFRSF1A, TNFRSF1B, IL1R1 (6), TBX21 (7) · negative: MUC1, EBF2, SLA2, MAP4K1,
  JMJD1C, GRAP, IKZF3
- **CRISPRa IL-2**: VAV1, CD28, LCP2, IL2RB, IL2 itself (5), EMP1, LHX6 ·
  negative: MAP4K1, LAT2, SLA2, GRAP, ARHGDIB, INPPL1, ITPKB
- **CRISPRi IFN-γ**: CD28, CD3G, VAV1, LCP2, ZAP70, CD3D, CD3E, CD247 depleted
  (all required), IFNG itself at 18 · enriched: CD5, CBLB, RNF40, MAP4K1, NFKB2
- **CRISPRi IL-2**: VAV1, LCP2, IL2 (5), CD3G, ZAP70, CD28, CD3D, CD3E depleted ·
  enriched: CBLB, CD5, EIF3D, RNF40, MAP4K1

Every directional claim the paper makes is reproduced. TBX21 is a positive hit in
IFN-γ (rank 7) and GATA3 a negative one (rank 35), and neither reaches
significance in IL-2 — the type-1/type-2 split the paper draws. The
MALT1–BCL10–TRAF6–MAP3K7–CHUK–IKBKG circuit is significant in CRISPRi IFN-γ
(ranks 43, 48, 9, 25, 59, 118) and reaches nothing comparable in the other three
screens, matching the paper's claim that this circuit is required rather than
limiting. The converse holds for the TNFRSF set: TNFRSF9 (rank 52) and IL1R1
(rank 6) are CRISPRa IFN-γ hits and sit at ranks 14,457 and 10,638 under CRISPRi
— the complementarity the paper's central argument rests on.

## Result 3 — head-to-head on recovery: level

| Screen | Spearman (effect) | BARCS hits | MAGeCK, published rule | MAGeCK, FDR only | shared |
|---|---|---|---|---|---|
| CRISPRa IL-2 | 0.85 | 157 | 477 | 792 | 143 |
| CRISPRa IFN-γ | 0.85 | 178 | 562 | 908 | 164 |
| CRISPRi IL-2 | 0.85 | 402 | 238 | 521 | 216 |
| CRISPRi IFN-γ | 0.84 | 329 | 211 | 413 | 190 |

Gene effects agree at Spearman 0.84–0.85 across all ~18,800 genes (0.91–0.94
restricted to either method's top 1000) with 82–84% sign agreement.

A 43-entry reference panel was fixed from the paper's text before either result
was looked at — each gene's biological role only, with the expected bin direction
derived from the modality, since CRISPRa and CRISPRi expect opposite shifts from
the same role. Recovery at matched list length:

| list length | BARCS | MAGeCK |
|---|---|---|
| top 50 | **36 / 43** | 34 / 43 |
| top 100 | 38 / 43 | **39 / 43** |
| top 250 | 39 / 43 | **40 / 43** |
| top 500 | 41 / 43 | **42 / 43** |
| each method's own rule | 40 / 43 | **41 / 43** |

BARCS leads at the very top of the list and MAGeCK by one gene from 100 onward.
A label-free precision proxy agrees that the two are level: counting KEGG
T-cell-receptor-signalling genes among each method's top 100, 250 and 500 across
the four screens, BARCS is ahead or tied in 6 of the 12 cells. **Neither method
demonstrably ranks this screen better than the other.**

BARCS misses three panel entries that MAGeCK calls — CBLB in both CRISPRa screens
(FDR 0.82 and 0.21) and SLA2 in CRISPRi IL-2 — and calls IKBKB in CRISPRi IFN-γ,
which MAGeCK misses.

**The hit counts are not a like-for-like comparison, and the reason matters.**
The published MAGeCK rule is a *conjunction*: FDR < 0.05 **and** median
|log2FC| > 0.5. The effect floor is not cosmetic. It removes 40%, 38%, 54% and
49% of the FDR-significant genes in the four screens, and on held-out controls it
halves the false positives, from 4/660 to 2/660. So the published rule does not
rest on its FDR being correctly calibrated — the effect floor is an empirical
guard that absorbs the anticonservatism. BARCS is read at a bare FDR < 0.05 with
no such floor, which is a different rule *shape*, not a different threshold. The
bare-FDR MAGeCK column above is the like-for-like statistical rule; against it
BARCS is the more conservative method in every screen.

## Result 4 — calibration, where the methods actually differ

Each screen's 165 control pseudo-genes are split into a calibration half and an
evaluation half by alternating pseudo-gene, and the two folds swap those roles.
Within a fold the calibration half is the only control set either method may tune
its null on — BARCS's control calibration, MAGeCK's `--control-sgrna` — and the
evaluation half is relabelled and scored alongside the ~18,800 real genes. Every
evaluation pseudo-gene is null by construction. Fitting a scale and evaluating it
on the same controls is not treated as validation.

Pooled over four screens and both folds, 660 held-out null pseudo-genes:

| Method | genes called | null genes at p < 0.05 | rate | matched expectation | ratio | null called at FDR < 0.05 |
|---|---|---|---|---|---|---|
| **BARCS** | 2,082 | 17 / 660 | 0.026 | 0.05 | **0.52** | **0** |
| BARCS, unmoderated | 1,677 | 24 / 660 | 0.036 | 0.05 | 0.73 | 1 |
| MAGeCK (published rule) | 2,810 | 82 / 660 | 0.124 | 0.10 | **1.24** | 2 |

The expectations differ because the statistics do: BARCS reports one two-sided
p-value, so 5% of null genes should fall below 0.05, whereas the published MAGeCK
rule reads the smaller of two near-complementary one-sided RRA tails, whose
first-order expectation is 10%. Against a strict 0.05 reference MAGeCK's 0.124
would read as 2.5× inflated; 1.24× is the fair figure, and the raw rate is in
`null_calibration.csv` for readers who prefer their own reference.

Two things follow. **Dispersion moderation is justified out-of-sample**: it adds
405 held-out discoveries (+24%) while cutting the null rate from 0.036 to 0.026
and the FDR-level false positives from 1 to 0. And **BARCS's shorter hit lists
are conservatism, not weakness** — at 0.52× nominal it is over-correcting by
about a factor of two, so the gap to MAGeCK's list length is a calibration
artefact with real power still recoverable, not a modelling deficit.

**Where the remaining conservatism sits.** `bb_gene_consistency()` sets its null
scale to `max(all-gene MAD, control-derived scale)`, so whichever term is larger
sets the scale for the whole screen. The control-derived term is a tail quantile
of at most 165 control-gene statistics and is the noisier of the two — across the
eight held-out fits it ranges 1.07 to 2.33 while the all-gene MAD, estimated on
~18,800 genes, is stable to three digits within a screen. A `max()` of a stable
estimator and a noisy one is dominated by the noisy one's upward excursions. A
gene-level null scale that does not hinge on a tail quantile of a few hundred
control genes is the concrete fix, and this screen is a ready-made test case for
one.

## What BARCS gives here that the published pipeline does not

These hold whether or not the hit lists differ, and they are the reason to run
the screen this way:

1. **No normalise-and-merge step.** Two separately sequenced pools are never put
   on a common scale; each guide is contrasted against its own libraries' totals.
2. **Donor is a model term**, not a pairing device, so it is estimated rather
   than differenced away and the design generalises past two donors.
3. **The unsorted bins are used.** The published test discards them; here they
   anchor the dispersion and raise the residual df from 3 to 5 before moderation.
4. **Both cytokine contrasts come from one fit per guide**, with a shared
   dispersion, rather than from four independent tests on four rebuilt tables.
5. **Raw counts throughout.** Nothing is rounded, rescaled, or renormalised
   between the GEO table and the likelihood.

## Limits — state these before any number is quoted

1. **Nontargeting controls are an optimistic null.** They carry no on-target
   effect *and* no locus perturbation, whereas a targeting guide against a
   phenotype-free gene still perturbs its locus. Every calibration number above
   is therefore a lower bound on the true null width, for both methods.
2. **One detail of the published recipe is inferred.** "Normalized to the total
   read count in each sample" is applied as MAGeCK's own total-count
   normalisation — count / library total × mean library total — which keeps the
   merged table on a read-count scale, the thing `--norm-method none` then
   expects. The merge itself needed no inference: the two library sets were
   checked and share no sgRNA id, the nontargeting sequences included.
3. **MAGeCK version.** 0.5.9.5 here; the paper used 0.5.9.2.
4. **The reference panel is small and literature-derived.** 43 gene × screen
   assertions taken from the paper's own text. It is a sanity check on ordering,
   not a truth set, and the one-gene differences it shows between methods are
   inside its resolution.
5. **The two hit rules are different rule shapes**, not two thresholds — see the
   end of Result 3. The held-out null check, which scores each method at its own
   rule, is the comparison that settles anything.
6. **Effects are logit-scale beta-binomial coefficients**, comparable to MAGeCK's
   median log2 fold change in sign and rank only, not in magnitude.
7. **A few genes are scored by MAGeCK and not by BARCS** — 125 in CRISPRa and
   74 in CRISPRi: `NO-TARGET` plus the genes left with fewer than three usable
   guides by `bb_gene_consistency()`. They are excluded from every joined
   comparison.
8. **Dispersion moderation is dataset-dependent.** It is applied here on the
   strength of the held-out check and the guide count; the Liu counter-example
   above is the reason it is not proposed as a package default.

## Reproducibility

```
Rscript scripts/prepare_schmidt_tcell.R              # GEO GSE174255 -> input/
Rscript examples/schmidt_tcell_barcs.R               # input/ -> output/
Rscript examples/schmidt_tcell_mageck.R              # input/ -> mageck/   (needs mageck)
Rscript examples/schmidt_tcell_null_calibration.R    # -> null_calibration.csv (needs mageck)
Rscript examples/schmidt_tcell_method_comparison.R   # -> comparison/, figure
```

`BARCS_NCORES` sets the guide-fitting parallelism (default 4); the eight
library-set fits take about 70 seconds in total at 8 cores. `MAGECK` sets the
path to the `mageck` executable if it is not on `PATH`. The 22 MB workbook, the
~18 MB of prepared inputs, and the guide-level tables all stay under git-ignored
paths; the compact tables under `comparison/` are committed. BARCS is
deterministic on this data.

## Suggested next steps

1. **Fix the gene-level null scale.** This screen is a clean test case: 165
   control pseudo-genes, two halves, and a null scale that moves 1.07–2.33
   between them. Any replacement — shrinking the control scale toward the
   all-gene MAD, pooling across related screens, or estimating it from a fitted
   two-component model rather than a tail quantile — can be scored here directly
   against the held-out protocol already in place.
2. **Decide the moderation default on the simulation grid**, not on real screens.
   `data/derived/barcs_null_calibration_grid_*` already sweeps guide count and
   overdispersion with known truth, which is what separates the Schmidt result
   (moderation helps) from the Liu one (it hurts).
3. **Use the supplementary CD4⁺ screens** (the paper's fig. S9: IL-2, IFN-γ and
   TNF-α in CD4⁺ T cells, one sequenced as two technical replicates). A
   technical-replicate arm would let the null be checked against something other
   than nontargeting controls, which is the limitation in Limits #1.

## Files

Input generator: `scripts/prepare_schmidt_tcell.R`. Shared design and helpers:
`R/schmidt_tcell.R`. BARCS fits: `examples/schmidt_tcell_barcs.R`. Published
pipeline: `examples/schmidt_tcell_mageck.R`. Held-out null check:
`examples/schmidt_tcell_null_calibration.R`. Comparison and figure:
`examples/schmidt_tcell_method_comparison.R`.

Committed under `results/schmidt_tcell/comparison/`: `method_concordance.csv`,
`null_calibration.csv`, `positive_control_panel.csv`, `tcr_pathway_precision.csv`,
and `<screen>_hits.csv` — every gene either method calls in that screen plus the
reference panel, with both methods' numbers side by side.

Regenerated, not committed: `input/`, `output/`, `mageck/`, `mageck_holdout/`,
`comparison/<screen>_gene_comparison.csv.gz`, and
`figures/schmidt_tcell_method_comparison.{png,pdf}`.
