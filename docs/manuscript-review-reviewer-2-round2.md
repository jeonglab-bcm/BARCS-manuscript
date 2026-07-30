# Referee report, second round — Reviewer 2

**Manuscript:** "BARCS: beta-binomial regression for multivariable CRISPR screen
designs" (revised; `865eeb9`)

**Recommendation:** Moderate revision. The revision is substantially and
honestly responsive, and the paper is now much closer to a claim its evidence
supports. The remaining work is mostly restoring material that the
restructuring dropped, plus one comparison that is not yet like-for-like.

All numbers below were checked against the manuscript source, `data/derived/`,
`results/`, `examples/`, and `R/bbreg.R`.

---

## Assessment of the revision

This is a genuinely responsive revision, and several changes go further than I
asked. Specifically:

- **Control calibration is now cross-fitted.** The IL2RA non-targeting rate is
  a deterministic five-fold held-out estimate (0.049), and the manuscript states
  outright why the in-sample version would be circular. `results/README.md` now
  records that policy. This was my first-round blocking point and it is fully
  resolved.
- **HT-29 is correctly demoted and the word "pseudoreplication" now appears in
  the paper.** Table S–HT29 includes all five deposited rows, and the text says
  plainly that "the full trajectory therefore provides no measured advantage
  over the endpoint on these two ranking metrics." That is the right conclusion
  and it takes some nerve to print it.
- **Sanson is now a model of how to report a null result.** Positive prevalence
  (63.4%) is stated, the gold-set realized FDP is reported (BARCS 0.0118 versus
  MAGeCK 0.0074 — i.e. the extra recall does come with more false calls, which
  is what I suspected), and the section concludes that neither model is
  established as better. Exactly right.
- **BARCS-MOD is now fully specified** (Equations 12–13: lowess prior scale,
  Smyth scaled-F prior df), including explicit disclosure that the `pmax`
  truncation is an implementation choice rather than a scientific constraint,
  and that all BARCS-MOD results are simulation ablations while the real-data
  benchmarks use unmoderated BARCS. That last sentence answers a question I
  had to infer from the source last time.
- **Prior art is cited** (`corncob`, `aod`, `VGAM`, `gamlss`) with a defensible
  positioning claim.
- **The Abstract and Discussion now state the losses**, including that
  alternative methods rank essential genes more strongly, and the CRISPulator
  section volunteers that edgeR-QL, DESeq2 and limma-voom "matched or exceeded
  BARCS ranking" in the 400-gene comparison.
- `docs/barcs-external-method-comparison.md` retracts its earlier over-reading
  in its own text ("When the general count models are included, F1 converges
  above a matched 0.005 and no method separates"). Correcting a supporting
  document rather than quietly leaving it is good practice.

**One correction to my previous report.** I wrote that the primary benchmark
tables were not deposited. That was wrong — I looked only in `results/`, and
`data/derived/` tracks the full Liang, CRISPulator and simCRISPR metric tables,
including per-seed values and standard errors. I withdraw that point, and it
changes several of my current asks from "run this" to "report the number you
already have."

I could not compile the manuscript (no TeX distribution in this environment),
so I have not verified the new arXiv-format build. The deleted `\repofile`
macro has no remaining uses and no `\includegraphics` target is missing, so
nothing obvious is broken.

---

## Major points remaining

### R1. The Liang calibration comparison is not like-for-like

This is now the paper's only head-to-head advantage claim on real data, so it
carries a lot of weight.

In `examples/liang_cas13_benchmark.R:364`, `bb_calibrate_controls()` is applied
to the BARCS guide statistics using the non-targeting controls. No equivalent
scaling is applied to edgeR-QL, DESeq2, limma-voom or MAGeCK-MLE. The reported
comparison — calibration error 0.025 for BARCS against 0.057–0.085 for the
others — therefore compares *calibrated BARCS* with *uncalibrated competitors*.

Control-based scaling of a *t* or Wald statistic is not specific to the
beta-binomial: the same 1,000 non-targeting guides could rescale any of the four
comparators. As it stands the result may show only that control calibration
works, which nobody disputes, rather than anything about the sampling model.

**Required.** Apply the same control-derived scaling to all five methods and
re-report the calibration errors. If the gap survives, the claim is much
stronger than it currently is. If it closes, say so — the paper is already
comfortable reporting that kind of outcome elsewhere, and the surviving claim
(a general-purpose depth-conditional coefficient model with an interpretable
calibration knob) would still stand.

### R2. Neither real-data design ablation now supports the central claim

The thesis is that fitting the extra design columns is worth doing. After this
revision, the two real-data tests of that thesis are:

1. **Liang three time points versus two.** The previous revision reported that
   adding day 7 raised average precision from 0.776 to 0.786 and FDR-0.10
   recall from 0.487 to 0.580. In the current version this ablation is gone
   entirely — I grepped for every one of those figures and for "two-endpoint"
   and found nothing. So the longitudinal section no longer contains any
   evidence that three time points beat two.

2. **IL2RA four bins versus outer bins.** Here the new "validation yield"
   column works against the argument. Four-bin BARCS yields 22/49 = 44.9%;
   outer-bin BARCS yields 18/34 = 52.9%. The marginal yield of the 15 extra
   discoveries the four-bin model makes is 4/15 = 26.7%, well below the
   outer-bin base rate. So by the paper's own newly introduced efficiency
   metric, adding the intermediate bins made the analysis *less* efficient. The
   text says only that four bins "recovered four regulators missed by its
   outer-bin fit," which is true but is the numerator without the denominator —
   and the very next sentence uses the yield metric to argue against MAUDE.

**Required.** Reinstate the Liang time-point ablation on the four
replicate-complete cell lines; the data are in `data/derived/` and the script
already supports it. And address the yield inversion in the IL2RA section
directly rather than around it. The honest framing is available and is not
weak: four bins recover more of the validated set (22 versus 18 of 26) at the
cost of lower precision per call, which is the ordinary sensitivity–specificity
trade and is a defensible thing to want. Stating it that way also removes the
appearance of applying the yield metric only where it flatters.

### R3. The generative-model calibration check is still not in the paper

`results/simulation_summary.csv` and the new `results/simulation_diagnostics.csv`
were both updated in this revision — clearly in response to the last round —
and still appear nowhere in the manuscript. They currently say:

| Statistic | BARCS *t* | BARCS moderated *t* | MAGeCK-MLE Wald |
|---|---|---|---|
| Gene null type I at 0.05 | 0.081 | 0.088 | 0.038 |
| Gene empirical FDR at nominal 0.05 | 0.130 | 0.091 | 0.049 |

with a guide-level null rate of 0.077 in `simulation_diagnostics.csv`. An
independent line of evidence points the same way: `control_scale_original` in
`data/derived/crispulator_facs_moi_10k_metrics.csv` is 1.25–1.32 across seeds,
so even in a clean four-replicate simulation the raw BARCS statistic needs
roughly 1.3× variance inflation to match its controls.

I want to be fair about two things. First, moderation improves the empirical FDR
substantially (0.130 → 0.091) — but the Supplement now states that the
real-data benchmarks use *unmoderated* BARCS, so the real-data analyses run in
the configuration with the worse number. Second, you tested my proposed
mechanism and it was wrong: `null_guides_at_rho_lower_boundary` is 0. I accept
that. But note that this simulation has 27 samples and (N − 1)ρ ≈ 60–180, so no
guide could plausibly sit at the boundary; the question is live in the actual
benchmark regime, where *m* is 8–12 and many guides are near-binomial.

**Required.** Report a nominal-versus-actual table in the Supplement across
*m* ∈ {4, 6, 8, 12} and a realistic range of ρ and guide abundance, at guide
and gene level, with and without moderation. Also report the fraction of guides
at the ρ̂ = 0 boundary in HT-29, IL2RA and Liang. A methods paper whose central
virtue is calibration should show its own null distribution, and a modest
anti-conservatism honestly reported and then repaired by moderation or control
scaling is a far better story than an unreported one.

### R4. The interaction simulation still has no valid comparator

Table 3's only external method returns F1 = 0.006, and the manuscript again
explains why the design is outside its scope. So the whole section remains a
BARCS-versus-BARCS denominator ablation, and the mechanism it identifies —
composition shift under widespread depletion — is not specific to the
beta-binomial.

edgeR-QL, DESeq2 and limma-voom are already run in the Liang benchmark and in
the 400-gene CRISPulator comparison, so fitting them here with an interaction
term is a small increment of work and would let the section separate what the
sampling model contributes from what the denominator contributes.

I would also restore the mechanistic detail this revision deleted: the median
estimated interaction of +0.18 among true zeros with an interquartile range
excluding zero, and the control inflation factors (1.99 and 2.24 for the
library denominator against 1.00–1.07 for the control denominators). Those
numbers were the evidence for the mechanism the section now asserts in one
clause.

### R5. Deleting `methods/1_model.tex` removed the paper's own central equation

The model subsection now gives Equation (1) and stops. Gone from the manuscript
entirely are:

- the hierarchy *K*<sub>i</sub> | *P*<sub>i</sub> ~ Binomial, *P*<sub>i</sub> ~ Beta;
- the definition ρ = 1/(κ + 1), so ρ appears in Equations (1), (10) and (11)
  without ever being defined as an intraclass correlation;
- **Var(*Y*<sub>i</sub>) = μ<sub>i</sub>(1 − μ<sub>i</sub>){(1 − ρ)/*N*<sub>i</sub> + ρ}**;
- the interpretation of the fitted coefficient (exp(β<sub>d</sub>) as an odds
  ratio per unit, and as an approximate log abundance ratio for rare guides).

The third of these is the mathematical statement of the paper's thesis. The
Introduction and Discussion both claim that the beta-binomial "separates
depth-dependent sequencing variation from heterogeneity among biological
libraries," and that separation *is* that equation. It should not be the one
display the manuscript omits. The fourth matters for a different reason: the
paper's primary output is a logit-scale coefficient, and it no longer tells a
reader how to interpret it.

**Required.** Restore all four to the Supplement's estimation subsection. This
is half a page and costs nothing.

### R6. Attribution regressions

Deleting `methods/0_cb2-to-barcs.tex` removed the paper's prior-art citations
along with its prose. The following are now in `main.bib` but cited nowhere:

- `williams1982extrabinomial` — the extra-binomial variance model that
  Equation (10)'s weight comes from;
- `baggerly2003differential`, `baggerly2004overdispersed` — the beta-binomial
  weighted test and the multi-group/continuous-covariate regression bridge;
- `li2014mageck` — MAGeCK itself. MAGeCK-MLE is the primary comparator in
  Liang, IL2RA, GSE70038, Sanson and HT-29, and the only MAGeCK citation left in
  the manuscript is a single `li2015mageckvispr` in the CRISPulator supplementary
  methods;
- `allen2019jacks`, `colic2019drugz` — JACKS and DrugZ are named in the audit
  section with no citation;
- `wang2019mageckflute` — the design-matrix rules the GSE70038 layout follows;
- `kolde2012rra`, `marioni2008rnaseq`, `spahn2017pinaplpy`,
  `eddelbuettel2014rcpparmadillo`.

There is an irony here: the revision added `corncob`, `aod`, `VGAM` and
`gamlss` in response to my last report while dropping Williams and Baggerly.
Please restore them, cite MAGeCK on first use in the Results, and cite JACKS
and DrugZ where they are named.

### R7. Comparator numbers were deleted rather than reported

The CRISPulator section now says the comparator set "was restricted to
MAGeCK-MLE and CRISPhieRmix; it therefore supports the moderation ablation, not
a universal method ranking." That caveat is correct, but the section reports no
comparator number at all, so a reader cannot see what the restriction conceals.
From `data/derived/crispulator_facs_moi_10k_metrics.csv`, three-seed means:

| Method | AUROC | Avg. precision | Realized FDP | F1 |
|---|---|---|---|---|
| BARCS-ST | 0.947 | 0.902 | 0.065 | 0.811 |
| BARCS-MOD | 0.954 | 0.919 | 0.066 | 0.847 |
| MAGeCK-MLE | 0.958 | 0.921 | 0.006 | 0.718 |
| CRISPhieRmix | 0.929 | 0.867 | 0.197 | 0.785 |

MAGeCK-MLE ranks better than BARCS-MOD on both AUROC and average precision
while calling an order of magnitude less liberally; CRISPhieRmix sits at 3×
the nominal error rate. The interesting content of this benchmark is that the
three methods occupy three different operating points, which is a more useful
message than the moderation ablation alone. The previous revision reported
these numbers and overclaimed on them; the fix is the caveat, not the deletion.

The same applies to the 400-gene comparison. The current sentence — the general
count models "matched or exceeded BARCS ranking but had higher realized FDP" —
is true at nominal thresholds, but `docs/barcs-external-method-comparison.md`
reports that at *matched* realized FDP no method separates above 0.005. Say that.

### R8. The Dempster audit needs a conflict statement and a narrower Abstract claim

Two separable issues.

**Disclosure.** CB<sup>2</sup> (`jeong2019beta`) is the corresponding author's
own prior method, and the audit rebuts a criticism of it. That is entirely
legitimate work, and the audit itself is careful and well documented. But the
Abstract now foregrounds the rebuttal, and readers will notice the relationship
whether or not the paper points it out. One sentence stating it — in the audit
subsection and in a competing-interests note — costs nothing and removes the
issue.

**Overstatement.** The Abstract says restoring full-library totals "eliminated
the reported CB<sup>2</sup> null discoveries" and concludes that the results
"support the beta-binomial framework in both CB<sup>2</sup> and BARCS." But the
Diagnosis column of the same table reports a largest full-total null *p* < 0.05
rate of 0.145 — nearly 3× nominal. Zero Benjamini–Hochberg discoveries at FDR
0.10 and a 0.145 nominal null rate are compatible, because BH responds to the
shape of the *p*-value distribution rather than its level. So the audit
establishes that the specific reported figure of 152 was an artifact of the
denominator contract; it does not establish that the null is calibrated, and it
cannot support "support the beta-binomial framework."

The Results subsection is appropriately careful ("These results do not
establish a winner against a complete ground truth"). Bring the Abstract into
line with it: the audit shows that pre-fit subsetting changed the effective
denominator and that the headline count does not survive correcting it, with
nominal-level calibration still unresolved.

A smaller symmetry point: the audit's argument is that upstream processing
changes the effective depth seen by a depth-conditional model. The Liang
analysis rests on ComBat-corrected, renormalized, rounded values. The revision
handles this much better than the last one — "sensitivity analysis" now appears
in the subsection title, the Discussion, and the Abstract's hedge — but one
clause acknowledging that the same caution applies to your own headline dataset
would close the loop.

---

## Minor points

1. **Report the deposited per-cell-line Liang table.** The macro-average hides
   strong heterogeneity that is scientifically interesting.
   `data/derived/liang_cas13_metrics_by_cell_line.csv` shows BARCS conservative
   in HEK293FT (null rate 0.029, average precision 0.662 against 0.728–0.753,
   FDR-0.10 recall 0.267 against 0.483–0.550) while the competitors run at
   0.13–0.20 in MDA-MB-231 and THP1. So the calibration advantage comes from two
   cell lines and the ranking deficit from a third. That is a more informative
   result than either average.
2. **Restore `essential_fdr_0_10_recall`.** It is deposited (BARCS 0.658 against
   0.800–0.813), it was reported in the previous revision, and it is the
   threshold an analyst actually uses. Omitting the metric with the largest gap
   while reporting two smaller ones invites exactly the reading the rest of the
   revision works to avoid.
3. **Add the deposited standard deviations to Table 3.**
   `data/derived/simcrispr_interaction_summary.csv` carries `sd` and `se` for
   every entry; with three seeds the table should show them.
4. **State how the three proxy-null trajectory examples were chosen**, and give
   the count of the reverse cases (proxy-null guides BARCS called and the others
   did not). Three hand-picked cases in an Abstract claim need their selection
   rule stated.
5. **Discuss proxy-null purity.** Unexpressed lncRNAs are a reasonable proxy,
   but TPM = 0 calls are noisy at low expression, and RfxCas13d has documented
   collateral activity, so the null set may contain guides with genuine fitness
   effects. That impurity rewards conservatism, which is the direction of the
   BARCS result. One or two sentences.
6. **Figure 1 axes.** A signed logarithmic effect axis combined with a vertical
   axis of log<sub>10</sub>{1 − log<sub>10</sub>(*p*)} is very hard to read; few
   readers will decode a double transform. Standard −log<sub>10</sub>(*p*) with
   clipping (as in the previous version, which marked clipped points with
   triangles) plus an inset for the dense region would serve better.
7. **Report the held-out IL2RA rate as a count and give the fold-wise scales.**
   The held-out 0.049 is 29/593, the same count as the all-controls fit, which
   is good evidence that the scale is stable across folds — worth saying
   explicitly. Please also report the five fold-wise scale estimates.
8. **Confirm the IL2RA fold structure.** Assigning folds cyclically by guide
   identifier is deterministic and reproducible, but if identifiers are ordered
   by plate or library position the folds will be structured rather than
   exchangeable. Please confirm, or add a permuted-assignment check.
9. **GSE70038 enrichment is close to a null result as reported.** "RNA
   metabolism" is the top pathway for both methods' exclusive sets in three of
   four coefficients, so the complementarity claim rests on GSC0131 alone, and
   the section concedes the Enrichr background problem. A decisive analysis is
   available at no cost: ask which method's exclusive hits are more enriched for
   the Hart reference essential genes you already use in the Sanson section.
   That gives an external criterion instead of a pathway-label comparison.
10. **The Methods promise Supplement content that does not exist.** The main
    Methods state that the Supplement "reports their complete preprocessing,
    design, and scoring definitions in the same order as the corresponding
    results," but there are no GSE70038 or Sanson methods subsections in the
    Supplement, and the restructuring dropped the GSE70038 design-matrix
    construction (the MAGeCKFlute Table 5 / Box 3 rules), its normalization, and
    the guide-to-gene aggregation used for Sanson. Add both subsections.
11. **The scale claim lost its evidence.** The Introduction's novelty claim is
    now that existing beta-binomial software does not "fit guide-level design
    matrices at screen scale," but the Software section that carried the runtime
    and convergence figures was deleted, so nothing in the manuscript supports
    it. Either restore a short scalability note (the runtime data are in
    `data/derived/crispulator_facs_multimethod_runtime.csv`) or drop the
    scale-based part of the claim. A direct timing comparison against `corncob`
    on 50,000 guides would settle it in one line.
12. **Losing the Limitations section costs more than it saves.** Several items
    that were material have no home in the current structure: linearity on the
    logit scale and depletion saturation (which Chronos models and a single
    slope does not); boundary fits for sparse or separated guides; the
    √*m*<sub>g</sub> independence assumption in the Stouffer combiner making
    gene statistics anti-conservative under shared-sequence correlation; and the
    plug-in nature of Equation (9), which does not propagate ρ̂ uncertainty. The
    Discussion's "two conditions" paragraph is good but covers only replication
    and calibration independence. A short Supplement limitations list would
    recover the rest.

---

## What would make this acceptable

1. Like-for-like Liang calibration — same control scaling for all five methods
   (R1).
2. Reinstate the Liang time-point ablation, and address the IL2RA yield
   inversion in the text (R2).
3. A nominal-versus-actual calibration table for the estimator itself, plus
   ρ̂-boundary fractions in the real benchmarks (R3).
4. Restore Equation Var(*Y*), the ρ definition, the beta-binomial hierarchy, and
   the coefficient interpretation (R5).
5. Restore the Williams, Baggerly and MAGeCK citations and cite JACKS/DrugZ
   (R6).
6. Report the CRISPulator comparator table with its caveat instead of omitting
   it (R7).
7. A competing-interests statement on the CB<sup>2</sup> audit and an Abstract
   claim matched to the 0.145 nominal null rate (R8).

Items 4–7 are restorations and reporting, not new analysis. Items 1–3 are real
work but bounded, and the scripts and deposited data already do most of it.
R4 (an interaction-design comparator) I would accept as a stated limitation if
the authors prefer, given that the section is now framed explicitly as an
ablation.

Closing note. Between the two rounds this manuscript has moved from claiming a
better-calibrated, higher-F1 method to claiming a genuine extension of the
beta-binomial screen model with an honestly characterized calibration–ranking
trade-off. That is a smaller claim and a much more defensible one, and the
sections rewritten most heavily — HT-29, Sanson, and the audit — are now among
the more scrupulous method comparisons I have read. The main risk in this
revision is the opposite of the last one: in tightening the claims, the
restructuring deleted supporting equations, citations and comparator numbers
that the paper needs in order to be read and checked. Putting those back is
mostly editorial, and would leave the manuscript in good shape.
