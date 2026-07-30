# Reviewer 2 — addendum on `hj/manuscript-minor-revision`

**Revision reviewed:** `abbbf39` on branch `hj/manuscript-minor-revision`
(1 commit ahead of `main` at `bb68802`, unmerged as of this note).

**Status:** Seven of the eight items in my round-3 acceptance list are done, two
of them more thoroughly than I asked. **Recommendation: accept with minor
revision**, subject to the three items below — one of which is a genuine
internal inconsistency created by the fix, and one of which is the most
important scientific finding in the paper and is currently confined to the
Supplement and Discussion.

For the record: `main` itself is unchanged since my round-3 report. The other
remote branches (`claude/main-text-figure-revision`, `manuscript`,
`hj/pixi-toolchain-and-readability`, `studio`, `hj/barcs-three-methods`,
`barc-manuscript`) contain no manuscript work newer than `main`.

---

## Round-3 items: what was done

**Verified** against the branch's own deposited data unless noted.

1. **simCRISPR per-seed restatement — done, and extended.** The section now
   reports all three seeds for both denominators and states plainly that "The
   three-seed means in Table 4 are not representative of a typical run because
   one library-denominator seed collapsed," with the corrected effect size
   "approximately 0.71 versus 0.91 rather than a near doubling." They also added
   the non-targeting per-seed values (0.798/0.924/0.889 and 0.899/0.933/0.920),
   which supports a claim I had not asked for and which is stronger than the
   original: the control denominator improved *every* seed. The Discussion and
   Abstract were updated to match.
2. **Null grid extended — done, and it settled the question against the
   authors' own earlier hypothesis.** The grid now spans
   (*N*−1)ρ ∈ {0, 5, 60, 180} at *N* = 100,000 — explicitly covering the
   continuous-dose simulation's 60–180 range — adds guides-per-gene ∈ {3, 5},
   adds a Gaussian-copula arm with within-gene correlation 0.4, and adds a
   split-control held-out gene column. All eight rows of the new Table S–grid
   **verified exactly** against `barcs_null_calibration_grid_summary.csv`. The
   finding is decisive and is reported as such: dispersion alone does *not*
   reproduce the 0.081 failure (largest independent-guide gene rate 0.074, and
   "error did not increase monotonically with dispersion"), whereas within-gene
   correlation drives gene-level error to 0.108–0.262 with guide-level error at
   or below nominal. This also retracts the previous Discussion claim that
   anti-conservatism grows with dispersion. Locating the failure in the
   combiner, and saying so, is the strongest single piece of work in this
   revision.
3. **One-guide pseudo-gene limitation — fixed rather than merely disclosed.**
   Non-targeting guides are now grouped into 156–159 pseudo-genes per cell line
   with guide counts drawn from the target-gene distribution
   (**verified**: `liang_cas13_control_pseudogene_sizes.csv`, sizes 2–10, mode 6),
   and the same assignments go to all five methods.
4. **CRISPulator aggregation — fixed.** Reverted to the MOI-0.20 three-seed
   means (0.902 → 0.919 AP, 0.811 → 0.847 F1, FDP 0.065/0.066), caption and
   Methods now consistently call 0.30 secondary.
5. **Cross-references and orphaned Methods — fixed.** **Verified** by static
   analysis: 22 `\ref`/`\eqref` calls, every figure and table now referenced
   from prose, no undefined reference, all three `\includegraphics` targets
   present, and **zero unused bibliography entries** (`dersimonian1986meta` and
   `spahn2017pinaplpy` removed with the sections that needed them). The
   BARCS-NORM / BARCS-RE / BARCS-EB / BARCS-GC specifications and the orphaned
   ten-scenario, 108-scenario and *R* = 1 CRISPulator Methods text are gone.
6. **Screen-scale runtime — partially addressed.** A number now exists: "Fitting
   the 50,000 guide-level regressions required 70.1–102.8 seconds per primary
   seed" (**verified** against the `guide_fit_seconds` column: 70.081, 102.787,
   92.55). A direct `corncob` comparison is still absent, so the Introduction's
   claim remains partly on assertion — but with a concrete screen-scale timing
   in the paper I no longer consider this blocking.
7. **Day-7 ablation — the concern largely dissolved under recalibration.**
   **Verified** per cell line: FDR-0.10 recall now improves in all four lines
   (0.733→0.783, 0.183→0.267, 0.617→0.717, 0.617→0.633). MDA-MB-231 remains
   marginally worse on AP (0.9441→0.9438) and on recall at 5% null FPR
   (0.967→0.950). The four-row table is still not shown, which is now a
   presentation preference rather than a substantive gap.

Not addressed: the IL2RA section title and the "strongest recovery-per-discovery
ratio among the four-bin models shown" sentence (round-3 point 8), and the
symmetric disagreement count (round-3 point 7 — the deposited
`liang_single_guide_disagreement_counts.csv` is unchanged, still only the two
asymmetric rows at *p* ≥ 0.20 and *p* ≥ 0.05).

---

## Three items before acceptance

### A1. The Abstract's calibration sentence is now stale and contradicts the revised Results

`sections/0_abstract.tex:9-12` still reads: "Applying the same
non-targeting-control scaling rule to BARCS and four competing count methods
produced **similar calibration for BARCS, edgeR-QL, DESeq2, and limma--voom**,
while the alternatives ranked essential genes more strongly."

That enumeration — four methods clustering, MAGeCK-MLE implicitly excluded — is
the *old* one-guide result, where MAGeCK sat at 0.0571 against 0.0296–0.0311.
Under the new aggregation-matched analysis the ordering **reversed**.
**Verified** against `liang_cas13_metrics_macro_average.csv` on this branch:

| Method | Mean abs. calibration error | AP | FDR-0.10 recall |
|---|---|---|---|
| MAGeCK-MLE | **0.0171** | 0.8771 | 0.7375 |
| limma-voom | 0.0178 | 0.8765 | 0.7667 |
| edgeR-QL | 0.0181 | 0.8745 | 0.7458 |
| DESeq2 | 0.0198 | 0.8739 | 0.7542 |
| **BARCS** | **0.0207** | **0.8382** | **0.6000** |

So BARCS is now last of five on calibration error, average precision, FDR-0.10
recall, and recall at matched 5% null FPR (0.867 versus 0.879–0.888). The
Discussion says this correctly — "with MAGeCK-MLE numerically closest to nominal
and BARCS lower in ranking and FDR recall" — but the Abstract's sentence
survives from the previous version and now implies the opposite grouping, and
the Results sentence ("aggregation-aware calibration removed the earlier
MAGeCK-MLE separation," `sections/results/2_liang-cas13.tex:39-41`) understates a
reversal as a removal.

**Fix.** Rewrite the Abstract sentence to match: all five methods calibrate to
0.017–0.021 once controls are aggregated to match the gene-level statistic, with
BARCS numerically last, and the alternatives also ranking essential genes more
strongly. In the Results, say "reversed" rather than "removed." This is the last
place in the manuscript where the front matter is ahead of the evidence, and
after three rounds of the authors correcting exactly this pattern it would be a
shame to leave it.

### A2. The correlated-guide finding is not carried into the benchmark sections it affects

This is the substantive point. The new grid shows that with within-gene
correlation *r* = 0.4, Stouffer gene-level type-I error reaches 0.232–0.262 for
five-guide genes and 0.124–0.144 for three-guide genes, while the guide-level
test stays at 0.027–0.048 (**verified**, all cells, against the deposited
summary). Aggregation-matched split-control scaling reduces this to 0.104–0.112
and 0.056–0.068 — still roughly 2× nominal in the five-guide case.

Every gene-level result in the manuscript uses that combiner: Liang, IL2RA,
CRISPulator, Sanson and GSE70038. Real guides targeting the same gene are
correlated through shared target biology, efficiency, and seed-based off-target
effects, so *r* = 0.4 is a plausible operating point rather than a worst case.
The implication is that the reported gene-level FDRs throughout the paper may be
anti-conservative by about a factor of two even after control scaling — and that
is the authors' own result.

**Verified** by grep: the finding appears only in `supp_calibration.tex` and one
Discussion paragraph. No benchmark section carries it, and
`sections/methods/3_inference.tex:127-131` still states the √*m_g* independence
assumption in the generic terms it used before the experiment existed.

**Fix.** Two sentences, not new analysis. (i) In
`sections/methods/3_inference.tex`, replace the generic caveat with the measured
one, citing Table S–grid. (ii) Add one sentence where gene-level FDRs are first
reported — the Liang section is the natural place — stating that the grid
quantifies the residual inflation under correlated guides and that gene-level
FDRs should be read with that factor in mind. This strengthens the paper: it
converts a hand-waved limitation into a measured one, and it is the natural setup
for the hierarchical gene model the Discussion already names as future work.

### A3. Two small items from round 3 that remain open

- `sections/results/3_facs-ordered-bin.tex:2-3` is still titled "A
  covariate-adjusted ordered phenotype **recovers validated regulators**" while
  the section shows the four-bin fit at *lower* validation yield than its own
  outer-bin fit (44.9% versus 52.9%) and MAUDE and Waterbear recovering more
  (25/26 and 24/26 versus 22/26). Suggest "An ordered-bin design trades
  precision per call for sensitivity."
- Line 75 of the same file still reads "BARCS therefore offered the strongest
  recovery-per-discovery ratio among the four-bin models shown," which is true
  only after the four-bin restriction; the MAGeCK outer-bin row in the same
  table is 60.0%. Add the outer-bin comparison to the sentence or drop the
  superlative.
- The symmetric disagreement count (BARCS *p* ≥ 0.05 with all three
  alternatives *p* < 0.05) is still not computed, so the 61-versus-1,242
  contrast in the Figure 1 caption still compares a *p* ≥ 0.20 criterion against
  a *p* ≥ 0.05 one. One extra row in
  `liang_single_guide_disagreement_counts.csv` and one clause in the caption.

---

## Assessment

This revision does something unusual: the authors ran the experiment I asked for,
it produced a result that removed their last remaining head-to-head advantage
and placed their method last on every metric of the Liang benchmark, and they
reported it — in the table, in the Results, and in the Discussion. The grid
extension likewise falsified their own stated mechanism (dispersion) and
identified a different one (the gene combiner) that is less flattering because it
implicates a component used in every benchmark in the paper.

What remains is A1 — a single stale Abstract sentence that now points the wrong
way — and A2, which is not a criticism so much as an instruction to give the
paper's best new finding the prominence it earned. Neither requires new
computation. With those two changes and the three cosmetic items in A3, I would
recommend acceptance.
