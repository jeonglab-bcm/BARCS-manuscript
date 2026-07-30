# Referee report, third round — Reviewer 2

**Manuscript:** "BARCS: beta-binomial regression for multivariable CRISPR screen
designs" (revised; `bb68802`)

**Recommendation:** Minor revision. The two substantive round-2 blockers are
genuinely fixed, the numerical fidelity of this version is high — I re-derived
essentially every table from the deposited data and found no mismatch — and the
claims are now calibrated to the evidence. What remains is one substantive
statistical gap, one result whose point estimate is an artifact of averaging,
and a set of editorial repairs the restructuring left behind.

Method: I checked each table against `data/derived/` and `results/`, read
`R/bbreg.R` and the analysis scripts, and ran static checks on the LaTeX and
bibliography. Numbers below marked **verified** are ones I recomputed or read
directly from a file.

---

## What this revision resolved

**From round 2:**

- **R1 (like-for-like Liang calibration) — fully fixed, and correctly.**
  `examples/liang_cas13_benchmark.R:168-193` defines one
  `calibrate_gene_controls()` (signed-*z*, 95th percentile of |control z| mapped
  to Φ⁻¹(0.975), `min_scale = 1`), and `add_method()` at line 476 applies it
  uniformly to all five methods before BH. Per-method scales are deposited to
  `data/derived/liang_cas13_control_scales.csv`. The result went against BARCS
  and is reported that way: calibration errors 0.0310 BARCS, 0.0296 DESeq2,
  0.0311 edgeR-QL, 0.0310 limma-voom, 0.0571 MAGeCK-MLE (**verified** exactly
  against `liang_cas13_metrics_macro_average.csv`), with the Discussion stating
  that "control scaling is therefore a general calibration device rather than
  evidence for one count likelihood." That is the right conclusion and it cost
  the paper its headline.
- **R2 (missing design ablations) — both fixed.** The day-7 ablation is back
  and deposited (`liang_cas13_timepoint_ablation.csv`): AP 0.8316 → 0.8382,
  FDR-0.10 recall 0.6083 → 0.6583, recall at 5% null FPR 0.8583 → 0.8667, null
  rate 0.0578 → 0.0709 (**verified**, all four). The IL2RA yield inversion is now
  stated head-on: "its overall yield was lower than the outer-bin fit (44.9%
  versus 52.9%)... they did not uniformly improve efficiency."
- **R3 (calibration of the estimator) — fixed and then some.** The new
  `sections/results/supp_calibration.tex` reports the failure openly ("exposes a
  finite-sample failure that a methods paper should not hide"), and the null grid
  in Table S–grid is **verified exactly**, all eight rows, against
  `data/derived/barcs_null_calibration_grid_summary.csv`. The real-data boundary
  audit is new: ρ̂ = 0 in 11/86,840 HT-29, 0/5,999 IL2RA, 1,738/280,541 Liang.
  That falsifies the boundary hypothesis I raised in round 1, on real data, and
  the paper says so.
- **R5 (deleted model equations) — restored.** The Binomial–Beta hierarchy,
  ρ = 1/(κ+1), Var(*Y*ᵢ) = μ(1−μ){(1−ρ)/*N*ᵢ + ρ} (Eq. 11), and the odds-ratio
  interpretation are all back in `sections/methods/2_estimation-irls.tex`.
- **R6 (attribution) — fixed.** Williams, Baggerly (2003, 2004), MAGeCK,
  JACKS, DrugZ, MAGeCKFlute, RRA, Enrichr, Smyth and RcppArmadillo are all now
  cited. **Verified:** exactly one bib entry is unused (`spahn2017pinaplpy`) and
  no cited key is missing from `main.bib`.
- **R7 (deleted comparator numbers) — fixed.** Table 3 now carries MAGeCK-MLE
  and CRISPhieRmix, the text states MAGeCK "slightly exceeded BARCS-MOD in AUROC
  and average precision," and — importantly — "After matching methods on realized
  FDP, no method separated from the others above an FDP of 0.005."
- **R8 (audit framing) — fixed.** A Competing interests statement is in
  `main.tex`, the audit subsection discloses the relationship in its second
  sentence, and both the Abstract and the audit's closing sentence now say
  nominal calibration is unresolved with the 0.145 rate attached.
- Minor items 1, 3, 4, 5, 7, 8, 9, 10, 12 from round 2 are all addressed
  (specialist comparators in the Abstract's framing, simCRISPR SDs, the
  disagreement-count selection rule, proxy-null purity, the 29/593 count and
  fold-wise scales, the permutation check, the Hart reference-essential
  follow-up, the GSE70038 and Sanson Methods subsections, and the Discussion's
  replication/calibration conditions).

**From round 1:** the in-sample calibration circularity, HT-29
pseudoreplication, Sanson prevalence and gold-set FDP, BARCS-MOD specification,
prior-art citation, and comparator-selection consistency are all now handled.

I have also verified the following are exactly right, having recomputed them:
Table S–Liang (all 20 rows), the Liang macro-averages, the IL2RA cross-fit
scales (1.3317–1.4176 and 1.3322–1.4168 permuted), the IL2RA 29/593 and 31/593
rates, all simCRISPR means and SDs, the simCRISPR true-zero medians (0.184, IQR
0.136–0.231 library; 0.004 non-targeting; 0.030 safe-harbor), and Table
S–calibration-diagnostic. All three `\includegraphics` targets exist and no
`\ref` is undefined.

---

## Major points remaining

### 1. The simCRISPR library-denominator result is an average of a bimodal outcome, and the reported effect size does not describe any run

**Claim at issue.** `sections/results/5_simcrispr-interaction.tex:60-63`:
"BARCS-MOD consequently reached F1 0.481 with the library denominator, compared
with 0.917 and 0.914"; Discussion: control denominators "nearly doubled F1."

**Verified evidence.** Per-seed F1 from
`data/derived/simcrispr_interaction_metrics.csv`, library denominator:

| Seed | BARCS-ST F1 | BARCS-MOD F1 |
|---|---|---|
| 20250724 | 0.039 | 0.011 |
| 20250725 | 0.709 | 0.680 |
| 20250726 | 0.716 | 0.754 |

Directional recall behaves the same way (0.020 / 0.549 / 0.559 for BARCS-ST).
The reported SDs — 0.389 and 0.409 on F1 — are **verified** against
`simcrispr_interaction_summary.csv` and are consistent with this. So the mean
0.481 is the average of one near-total collapse and two runs at ~0.70–0.75; no
simulation produced anything close to 0.481. On two of three seeds the honest
contrast is roughly 0.71 → 0.91, not 0.48 → 0.92.

This is a case where adding the SDs I asked for in round 2 exposed something the
means concealed, so the fix is straightforward but it does change the claim.

**Fix.** Report the three per-seed values in the table or its caption; state the
effect as "F1 0.68–0.75 on two seeds and 0.01–0.04 on the third, against
0.91–0.92 under either control denominator"; diagnose what makes seed 20250724
collapse (the natural candidate is the control scale, given the 1.99/2.24 versus
1.00–1.07 figures already reported — check whether that seed's scale is extreme);
and drop "nearly doubled" from the Discussion in favour of the range. Three
seeds is too few for a mean when the outcome is this dispersed — either run more
seeds or present the distribution.

### 2. The null grid does not cover the dispersion regime in which the observed failure occurs, and the Discussion extrapolates past it

**Claim at issue.** `sections/6_discussion.tex:35-38`: "The null grid further
showed that unmoderated inference is conservative in the smallest near-binomial
designs and can become modestly anti-conservative as sample size and dispersion
increase."

**Verified evidence.** From `barcs_null_calibration_grid_summary.csv`, the grid
uses *N* = 10⁶ with ρ ∈ {0, 5×10⁻⁷, 5×10⁻⁶}, so (*N*−1)ρ ≤ 5. From
`examples/simulation.R:46-47`, the continuous-dose simulation that yields
type-I 0.081 uses ρ = 0.0015 with library sizes 40,000–120,000, so
(*N*−1)ρ ≈ 60–180 — roughly 30× beyond the grid's largest value. Within the grid
the unmoderated statistic is *conservative* everywhere except mildly at *m* = 12
(gene type-I max 0.062; guide max 0.057; at *m* = 4 it is 0.000–0.017 guide and
0.000–0.026 gene). The grid therefore does not contain a cell that reproduces
0.081, and the two results are never reconciled: the Supplement reports both on
the same page and the Discussion asserts a trend in ρ that the grid cannot test.

The paper is honest that the grid is "deliberately near-binomial" — the
`rho_zero_fraction` column runs 0.00–0.64, against 0.013%–0.62% in the three
real fits — but that candour is exactly the problem: the grid probes the regime
where boundary truncation dominates, not the regime real screens occupy.

**Fix.** Extend the grid to ρ ≈ 10⁻⁴ and 10⁻³ at the same *m* values, which is
the regime of the dose simulation and of the real fits. That will either
reproduce the 0.081 (in which case you have located it and can say so) or show
it is specific to the dose simulation's other features. Either outcome is
publishable; the current pair of unreconciled results is not.

### 3. The gene combiner is the untested link between the two calibration results

**Claim at issue.** `sections/methods/3_inference.tex` states that the Stouffer
step "give[s] equal weight to valid guides and use[s] the independence reference
√*m_g*; shared sequence artifacts or other correlation can therefore make the
gene statistic anti-conservative." Nothing tests this.

**Verified evidence.** In the grid, gene type-I exceeds guide type-I in most
cells (unmoderated *m* = 4: guide 0.000–0.017, gene 0.000–0.026; *m* = 12: guide
0.014–0.057, gene 0.018–0.062), so the combiner is already inflating relative to
the guide test even with independent guides. The dose simulation, where type-I is
0.081, uses 3 guides per gene (`examples/simulation.R:15-16`) against the grid's
5, and all guides within a gene there share a common true effect and a common
baseline draw. The combiner is therefore the most likely single explanation for
the discrepancy in point 2, and it is the one factor the grid holds fixed.

**Fix.** Add guides-per-gene as a grid factor and add one arm with
within-gene correlated guide effects. Report guide-level and gene-level type-I
separately, which the deposited grid already does — the columns exist, they just
are not used to attribute the failure.

### 4. The like-for-like calibration is estimated on one-guide statistics and applied to multi-guide statistics

**Claim at issue.** `sections/methods/7_benchmark-liang.tex:28-30`: "every
non-targeting guide was treated as a one-guide pseudo-gene. After each method's
native gene summary, the 95th percentile of the absolute signed-normal control
statistic was mapped to the two-sided standard-normal 0.05 cutoff."

This is the right way to make the five methods comparable, and I want to be
clear it is a genuine improvement. But a scale estimated from *m_g* = 1
statistics cannot correct miscalibration *introduced by* aggregating *m_g* > 1
correlated guides — which, by point 3, is where a substantial part of the
miscalibration lives. The essential-gene positives are multi-guide genes; the
controls calibrating them are not. The consequence is that the residual
calibration errors in Table S–Liang (0.0296–0.0571) measure the guide-level fit
plus whatever the combiner adds, corrected only for the former.

**Fix.** Either build the pseudo-genes from random *m_g*-guide bundles of
non-targeting controls matched to the gene-level *m_g* distribution, or state
this limitation explicitly where the rule is defined. The second option is
acceptable; the current silence is not, because the whole point of the rule is
that it makes the methods comparable.

### 5. The CRISPulator headline pools the primary and supplementary designs, contradicting the Methods

**Claim at issue.** `sections/results/4_crispulator.tex:6-7` says "Three
prespecified seeds"; the Table 3 caption says "three-seed mean across simulated
multiplicities of infection"; `sections/methods/5_benchmark-crispulator.tex:17-18`
designates MOI 0.20 "the main analysis" and 0.30 "the supplementary analysis."

**Verified evidence.** Recomputed from
`data/derived/crispulator_facs_moi_10k_metrics.csv`, low-bulk-high design:

| Aggregation | BARCS-ST AP / FDP / F1 | BARCS-MOD AP / FDP / F1 |
|---|---|---|
| MOI 0.20 only (n=3) | 0.9023 / 0.0651 / 0.8107 | 0.9188 / 0.0661 / 0.8468 |
| Pooled both MOI (n=6) | 0.9026 / 0.0698 / 0.8171 | 0.9195 / 0.0680 / 0.8487 |
| **Manuscript** | **0.903 / 0.070 / 0.817** | **0.919 / 0.068 / 0.849** |

So the reported values are the pooled n = 6 means, not three-seed means at the
main MOI. This also explains, silently, why the numbers moved from the previous
revision (which reported 0.902 / 0.065 / 0.811 and 0.919 / 0.066 / 0.847 — the
MOI-0.20 values). Nothing is wrong with pooling, but the text says "three
seeds," the Methods call 0.30 supplementary, and the change between revisions is
unflagged.

**Fix.** State n = 6 in the caption and the text, or report MOI 0.20 as primary
with 0.30 as a sensitivity row, matching the Methods. Either way, note that the
aggregation changed from the previous version.

---

## Moderate points

### 6. The day-7 ablation is mixed across cell lines and only the macro-average is shown

**Verified** from `data/derived/liang_cas13_timepoint_ablation.csv` (days 0+14 →
days 0+7+14):

| Cell line | AP | Recall @5% null FPR | FDR-0.10 recall |
|---|---|---|---|
| HAP1 | 0.8754 → 0.8854 | 0.883 → 0.900 | 0.733 → 0.783 |
| HEK293FT | 0.6616 → 0.6659 | 0.700 → 0.700 | 0.183 → 0.267 |
| MDA-MB-231 | 0.9441 → **0.9438** | 0.967 → **0.950** | 0.800 → 0.867 |
| THP1 | 0.8452 → 0.8576 | 0.883 → 0.917 | 0.717 → **0.717** |

The gain is carried by HAP1 and HEK293FT; MDA-MB-231 is slightly *worse* on two
of three metrics and THP1's FDR recall is unchanged. With four units, mixed
directions and no paired test, "improved essential-gene recovery modestly" is
about as strong as the data allow — but the per-cell-line breakdown is deposited
and should be shown rather than summarized away, especially since the paper
already shows the per-cell-line table for the method comparison.

**Fix.** Add the four-row ablation table to the Supplement and either a paired
test across cell lines or one sentence naming the two lines that improve and the
one that does not.

### 7. The 1,242-versus-61 asymmetry uses asymmetric thresholds and lives only in a caption

**Verified** from `data/derived/liang_single_guide_disagreement_counts.csv`: the
61 count uses "BARCS *p* ≥ 0.20; each alternative *p* < 0.05" and the 1,242
count uses "BARCS *p* < 0.05; each alternative *p* ≥ 0.05." The thresholds
differ, so the 20-fold ratio is not a like-for-like comparison, and the
symmetric count (BARCS *p* ≥ 0.05 with all three alternatives *p* < 0.05) is not
reported.

Credit where due: the caption reports both directions and says explicitly that
the panels "illustrate the trajectory structure of one disagreement direction,
not its prevalence." But the direction that runs against BARCS — 1,242
proxy-null guides that only BARCS calls — appears nowhere in the main text, and
it is arguably the most informative single number in the section.

**Fix.** Report the symmetric count, put both numbers in the text, and add one
sentence on what a 1,242-guide excess implies at guide level given that the
gene-level calibration errors are near-identical across methods. That
reconciliation is interesting in its own right.

### 8. Section titles and one comparative sentence overstate relative to their own sections

- `sections/results/3_facs-ordered-bin.tex:2-3`, "A covariate-adjusted ordered
  phenotype **recovers validated regulators**" — the section's own analysis shows
  the four-bin fit at lower validation yield than its outer-bin fit (44.9% vs
  52.9%) and MAUDE and Waterbear recovering more (25/26 and 24/26 vs 22/26).
  Suggest something like "An ordered-bin design trades precision per call for
  sensitivity."
- Same file, "BARCS therefore offered the strongest recovery-per-discovery ratio
  among the four-bin models shown" — true only after the four-bin restriction;
  the MAGeCK outer-bin row *in the same table* is 60.0% (**verified** against
  `results/waterbear_facs/benchmark_metrics.csv`: 30 discoveries, 18/26).
  Restricting the comparison class to win is the pattern this revision otherwise
  avoids. State the outer-bin numbers in the same sentence.

### 9. The Discussion generalizes a calibration result that holds in half the cell lines

`sections/6_discussion.tex:15`: "MAGeCK-MLE remained farther from nominal on the
proxy-null set." **Verified** per-cell-line calibration error, MAGeCK vs BARCS:
HAP1 0.0780 vs 0.0188 and MDA-MB-231 0.0901 vs 0.0363 (BARCS better), but
HEK293FT 0.0151 vs 0.0203 and THP1 0.0450 vs 0.0487 (MAGeCK better). The macro
gap comes entirely from two lines. The Results paragraph does point at HAP1 and
MDA-MB-231; the Discussion should carry the same qualifier — "in two of four
cell lines."

### 10. The Introduction's sole novelty claim has no support anywhere in the paper

`sections/1_introduction.tex:22-28` positions BARCS against `corncob`, `aod`,
`VGAM` and `gamlss` on the grounds that they "do not provide a CRISPR-oriented
workflow that preserves immutable full-library totals, fits guide-level design
matrices, and connects the result to established screen outputs." Since the
Lemma and the equivalence Proposition were removed in an earlier revision, and
the Software section with the runtime and convergence figures was removed in the
revision before this one, this sentence is now the entire novelty claim — and
nothing in the manuscript substantiates it. The only surviving trace is one
sentence about RcppArmadillo at `sections/methods/2_estimation-irls.tex:71-74`.

The evidence exists: `data/derived/crispulator_facs_multimethod_runtime.csv` and
`crispulator_facs_multimethod_runtime_summary.csv` are committed.

**Fix.** Add a short Supplement paragraph with (a) a timing comparison against
`corncob` fitted with fixed totals on the 50,000-guide CRISPulator library, and
(b) whether `corncob` with the same design and totals reproduces the BARCS
coefficients. Two numbers and one concordance statistic would settle it. If you
prefer not to run it, drop "at screen scale" and frame the contribution purely
as the screen-specific contract.

---

## Minor points

11. **Cross-referencing is broken, and Figure 1 lost its only citation this
    round.** **Verified** by static analysis: the manuscript contains 16
    `\ref`/`\eqref` calls in total, and the following labels are never
    referenced from any prose: `fig:liang`, `tab:liang-cell-lines`,
    `tab:simcrispr`, `fig:simcrispr`, `tab:chronos`, `tab:sanson`,
    `tab:chronosaudit`, `eq:barcs-model`. `sections/results/2_liang-cas13.tex`
    had one `\ref{fig:liang}` at `865eeb9` and has none now. So two of three
    figures and four tables — including the main Liang table and the paper's
    defining equation — are never pointed to. Easily fixed and will be flagged
    immediately at submission.
12. **About 130 lines of orphaned Methods.** BARCS-NORM, BARCS-RE
    (`BARCS-partial`), BARCS-EB and BARCS-GC are fully specified in
    `sections/methods/3_inference.tex` (Eqs. `eq:genenormal`, `eq:genenormalt`,
    `eq:generandomeffects`, `eq:genere`, `eq:geneeb`, `eq:guideconsistency`,
    `eq:gcmad`, `eq:gccontrol`, `eq:gcz`) and **verified** to appear in no
    Results section. Likewise
    `sections/methods/5_benchmark-crispulator.tex:45-68` still specifies a
    ten-scenario sensitivity analysis (50 simulations), an available
    108-scenario factorial, and the *R* = 1 / BARCS-GC boundary design, none of
    which is reported. Either report them — `README.md` already carries the
    numbers (AP 0.917 / 0.777 / 0.842 / 0.903 for original / NORM / partial /
    EB) — or move the specifications to a software appendix and say where the
    results live.
13. **CB² is still not a comparator in any benchmark.** **Verified:** in
    `sections/results/7_sanson-a375.tex` and
    `sections/methods/9_benchmark-sanson.tex`, CB² appears only as the source of
    the bundled Sanson data. The paper's framing claim is that BARCS
    "generalizes" CB² and "retains this sampling principle"
    (`sections/2_results-intro.tex:1-6`), and since the equivalence Lemma was
    removed, nothing checkable supports it. Running
    `bbreg(~ group)` against `CB2::measure_sgrna_stats()` on the Sanson guides —
    a screen already in the repository, with CB² already an `Imports` dependency
    in `DESCRIPTION` — and reporting the discrepancy distribution would close
    this in one paragraph.
14. `spahn2017pinaplpy` is the one unused bib entry (**verified**).
15. `DESCRIPTION` lists `CB2` under `Imports`, so BARCS depends on CB² at
    install time. That relationship disappeared from the manuscript when the
    Software section was deleted and is worth one sentence.
16. The Table 3 caption ("Each method is evaluated at its nominal gene-level
    threshold; the resulting operating points are therefore not matched on
    realized FDP") is good practice and should be echoed in Table S–Liang and
    Table 2, where the same caveat applies.

---

## What no one has raised yet

17. **The moderated statistic is referred to *t*_{ν+d₀} with d₀ estimated from
    the same screen.** `sections/methods/3_inference.tex` gives
    ν_mod = ν + d₀ with d₀ from the Smyth scaled-*F* moment estimator, and
    reports d₀ ≈ 42.7–45.8 in CRISPulator. Adding ~43 degrees of freedom
    estimated from the data to a residual ν of 7 is a large borrowing, and the
    reference distribution treats d₀ as known. This is standard limma practice,
    so I am not calling it wrong — but limma's *t* reference is validated for
    per-gene residual variances, whereas here the moderated quantity is a
    beta-binomial *dispersion inflation*, and the grid (point 2) shows moderation
    at *m* = 4 moving the gene rate from 0.000–0.026 up to 0.034–0.064, i.e.
    over-shooting in the worst cell. One sentence acknowledging that d₀ is
    plugged in, and a note that the grid is the evidence for how much that
    matters, would cover it.
18. **The IL2RA discovery counts come from the all-controls scale while only
    the NTC rate is cross-fitted.** The paper says this
    (`sections/results/3_facs-ordered-bin.tex:26-28`: "The production gene
    analysis used all controls to estimate one scale"), so it is disclosed, not
    hidden. But the fold-specific scales span 1.332–1.418 (**verified**), a 6%
    range, and no one has checked whether the 49 discoveries and 22/26 recovery
    are stable across that range. Recomputing the table at the minimum and
    maximum fold scale is three lines of code and would either confirm
    robustness or reveal that the headline count is scale-sensitive.
19. **`min_scale = 1` makes every calibration in the paper one-sided.** This is
    now applied uniformly across methods in Liang, so the comparison is fair —
    but it means no reported "calibration error" can ever reflect a correction
    for conservatism. In HEK293FT, BARCS, edgeR-QL and DESeq2 all sit at
    0.025–0.030 (**verified**), i.e. conservative, and the truncation leaves them
    there. The metric therefore measures "distance from nominal after an
    anti-conservatism-only correction," which is worth naming once where the rule
    is defined.

---

## What would make this acceptable

1. Restate the simCRISPR denominator effect from the per-seed values (point 1).
2. Extend the null grid to ρ ≈ 10⁻⁴–10⁻³ and add guides-per-gene as a factor, so
   the 0.081 is either reproduced or localized (points 2 and 3).
3. State the one-guide-pseudo-gene limitation of the like-for-like calibration
   (point 4).
4. Fix the CRISPulator aggregation description to n = 6, or report MOI 0.20 as
   primary (point 5).
5. Show the per-cell-line day-7 ablation and the symmetric disagreement count
   (points 6 and 7).
6. Retitle the IL2RA section and qualify the two comparative sentences (points 8
   and 9).
7. Support or drop the "at screen scale" novelty claim (point 10).
8. Repair the cross-references and resolve the orphaned Methods (points 11
   and 12).

Items 3–8 are editorial or one-paragraph additions. Item 1 is a restatement from
data already in the repository. Item 2 is the only real computation, and it is a
rerun of a script that already exists with two extra ρ values.

---

## Closing

Across three rounds this manuscript has moved from claiming a better-calibrated,
higher-F1 method, through a calibration–ranking trade-off, to what it now
claims: a multivariable extension of the beta-binomial screen model, with the
calibration advantage explicitly retracted once the comparison was made fair,
and a finite-sample failure reported in a table rather than omitted. The
Abstract now ends "not evidence of uniform superiority over negative-binomial
methods." Very few method papers are willing to write that sentence, and the
paper is stronger and more useful for it.

The remaining substantive gap is narrow and specific: the calibration study
probes the wrong dispersion regime to explain the one failure the paper itself
reports, and the gene combiner — the most likely culprit — is the one factor the
grid holds fixed. Closing that would make the Supplement's calibration section
genuinely conclusive instead of merely honest. Everything else on my list is
repair work from the restructuring: cross-references, orphaned Methods, one
aggregation label, and one unsupported sentence in the Introduction.

*Note on process:* I attempted to parallelize this review across six specialist
readers with adversarial verification of each finding. The subagent environment
failed on its first tool call and the agents returned nothing; the synthesizing
agent correctly refused to produce a report from no evidence rather than
inventing citations. Every finding above is therefore from my own direct
verification against the files named, and I have marked which numbers I
recomputed. The one number I initially cited without recomputing — the simCRISPR
true-zero median and IQR — I subsequently confirmed as 0.184 (0.136–0.231).
