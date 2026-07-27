# BARCS compared with other methods

This report separates two questions that can otherwise give contradictory
answers:

1. **Ranking:** are truly active genes near the top?
2. **Decision calibration:** when a method reports FDR below 0.10, is the
   realized false-discovery proportion actually near or below 0.10?

No method wins both questions in every benchmark.

**Scope note.** The manuscript compares BARCS against CRISPR-specific gene
callers only — MAGeCK-MLE and CRISPhieRmix — because that is the choice a
screen analyst actually faces. The general RNA-seq count models (edgeR-QL,
DESeq2, limma-voom) were dropped from it. Their results are kept here as the
supporting record, since they are what established that a fixed nominal cutoff
means very different things to different methods, and that finding still
motivates the calibration analysis the manuscript reports. DESeq2 also remains
in the genome-scale pipeline as CRISPhieRmix's documented input, supplying
guide-level log2 fold changes rather than competing as a gene caller.

## CRISPulator FACS simulations

The comparison draws on the same primary low + bulk + high simulations as the
four-method BARCS benchmark: ten one-at-a-time parameter scenarios and five
fixed seeds. Every method receives the same guide counts, ordered-phenotype
score, replicate adjustment, and gene truth.

**Analysis scope.** The headline comparison covers the 45 runs from the nine
settings with more than one independent screen replicate. The one-replicate
setting is a prespecified diagnostic boundary, not a supported operating
point, and it is reported separately below. At $R=1$ the low + bulk + high
design leaves a single residual degree of freedom per guide, so BARCS-original
makes no calls in any seed; pooling those five zeros into the mean reports a
boundary failure as if it were average performance. The scope is set by
`min_replicates` in
`examples/crispulator_facs_external_head_to_head_aggregate.R`, and every
derived table carries a `scope` column.

| Method | Average precision | AUROC | Directional recall at FDR 0.10 | Realized FDP | F1 | Negative-control $p<0.05$ |
|---|---:|---:|---:|---:|---:|---:|
| BARCS-original | 0.880 | 0.932 | 0.709 | 0.094 | 0.780 | 0.060 |
| BARCS-NORM | 0.727 | 0.880 | 0.175 | 0.065 | 0.272 | **0.048** |
| BARCS-partial | 0.803 | 0.899 | 0.419 | 0.047 | 0.559 | 0.014 |
| BARCS-EB | 0.864 | 0.923 | 0.502 | **0.014** | 0.643 | 0.013 |
| MAGeCK-MLE | 0.883 | 0.937 | 0.651 | 0.064 | 0.756 | 0.039 |
| edgeR-QL | 0.897 | **0.941** | 0.835 | 0.250 | 0.787 | 0.128 |
| DESeq2 | **0.897** | 0.940 | **0.837** | 0.243 | 0.793 | 0.122 |
| limma-voom | 0.897 | 0.941 | 0.828 | 0.229 | **0.794** | 0.119 |

Bold identifies the largest ranking, recall, or F1 value, the smallest
realized FDP, and the negative-control rate closest to its nominal 0.05
target. A very small FDP is conservative rather than automatically optimal;
it should be read together with recall.

The result is a real trade-off:

- edgeR-QL, DESeq2, and limma-voom rank simulated active genes slightly
  better and recover many more at their nominal threshold.
- Their nominal FDR 0.10 calls are anti-conservative here: mean realized FDP
  is 0.23--0.25, and 0.12--0.13 of negative controls have \(p<0.05\).
- BARCS-original's ranking deficit is small but real: paired average
  precision is 0.017 below edgeR-QL (95% interval 0.013 to 0.021). Its
  realized FDP is 0.156 lower (0.135 to 0.176).
- **On F1 the three count models are not distinguishable from
  BARCS-original once the boundary setting is removed.** Paired F1
  differences are $+0.008$ for edgeR-QL (95% interval $-0.015$ to $+0.030$),
  $+0.013$ for DESeq2 ($-0.010$ to $+0.035$), and $+0.014$ for limma-voom
  ($-0.008$ to $+0.037$). Pooling the one-replicate setting inflated the same
  three differences to $+0.063$, $+0.067$, and $+0.070$, all nominally
  significant. The apparent F1 gap was an artifact of the boundary.
- BARCS-NORM is well calibrated but underpowered. Estimating one standard
  deviation from roughly five guide beta values gives only about four
  reference degrees of freedom.
- BARCS-EB is the safest method, but the cost is substantial recall.
- MAGeCK-MLE is well calibrated and indistinguishable from BARCS-original in
  ranking. On paired runs, its average-precision difference from
  BARCS-original is $+0.0022$, with a 95% interval from $-0.0042$ to
  $+0.0086$; its F1 is lower by 0.0241 (95% interval $-0.0404$ to
  $-0.0079$).

### Why the ranking gap persists

The gene combiner is not the source of the difference. edgeR-QL, DESeq2, and
limma-voom are fitted guide by guide and then aggregated with the same
signed-$z$ rule as BARCS-original, so the comparison isolates the guide-level
count model.

Two properties of that model explain the two halves of the trade-off:

- `.bb_estimate_rho()` in `R/bbreg.R` fits the beta-binomial overdispersion
  separately for each guide, from that guide's own residuals. edgeR-QL,
  limma-voom, and DESeq2 shrink each guide's dispersion toward a trend fitted
  across the whole library. With $2R-1$ residual degrees of freedom per guide
  --- seven at the four-replicate baseline --- the unshared estimate is noisy,
  and that noise propagates into the standard error and degrades the ranking.
  This is the source of the residual average-precision gap.
- Shared shrinkage assumes a common dispersion--abundance trend. CRISPulator
  generates heterogeneous guide noise, so genuinely noisy guides are shrunk
  toward the pool and their variance is understated. That is what produces
  the 0.23--0.25 realized FDP and the 0.12--0.13 negative-control rate.
- `bb_calibrate_controls()` is a one-way ratchet: its scale is
  `max(min_scale, empirical/reference)` with `min_scale = 1`, so it can only
  inflate standard errors. It buys calibration at a cost in power. Because
  every guide shares the same residual degrees of freedom, this single global
  scale is rank-preserving at guide level and therefore does *not* contribute
  to the average-precision gap. Note also that the external methods receive
  no comparable step; part of the FDP difference is this method component
  rather than an intrinsic property of the count models.

### F1 across nominal FDR thresholds

For a called gene set, the standard binary F1 score is

$$
\mathrm{F1}=\frac{2TP}{2TP+FP+FN}.
$$

Here a true positive is any simulated active gene passing the nominated
gene-FDR threshold; direction is reported separately. The threshold scan uses
0.10, 0.05, 0.01, 0.005, and 0.001 in each of the 45 multi-replicate runs.

The scan separates calibration from power. BARCS-original's realized FDP
tracks its nominal threshold at every point (0.000, 0.004, 0.009, 0.054,
0.094 against nominal 0.001 through 0.10). The three count models overshoot
at every point, by 2.3- to 2.5-fold at nominal 0.10 and by 10- to 14-fold at
nominal 0.001. So the nominal FDR knob means what it says for BARCS and
does not for them.

Matching on *realized* FDP rather than the nominal label, the power gap is
smaller but does not vanish. At realized FDP near 0.06, BARCS-original scores
F1 0.753 (nominal 0.05) against 0.817 for edgeR-QL and 0.819 for DESeq2 (both
nominal 0.01). Thus the general count models retain a real power advantage
after moving to a threshold with comparable empirical FDP in this simulator;
their advantage is not explained entirely by the anti-conservative 0.10
operating point.

This is a diagnostic, not permission to select a method-specific threshold
from the known truth in new data. The selected 0.01 threshold would need
independent calibration or negative controls before prospective use.

At the four-replicate baseline, edgeR-QL has the highest mean average
precision (0.932), whereas BARCS-original has 0.917. limma-voom has the
highest F1 (0.832) and BARCS-original has 0.828, but their realized FDPs are
0.218 and 0.086, respectively.

### What the count models catch that BARCS misses

Aggregate metrics do not say whether the count models find *different* genes or
merely the *same* genes at a looser threshold. `examples/crispulator_facs_miss_cases.R`
answers that on the committed baseline realization in
`data/derived/crispulator_facs/` (seed 20250724, MOI 0.25, 90% high-quality
guides, 400 genes, four replicates; the MD5s match the head-to-head
provenance). BARCS is refit from scratch and reproduces the committed
head-to-head row exactly --- 59 calls, AUROC 0.9614, average precision 0.9225,
realized FDP 0.0678. edgeR-QL, DESeq2, and limma-voom are refit here rather
than reused, and land within one to two calls of their stored values.

**BARCS calls a strict subset.** Zero genes are called by BARCS and missed by
edgeR-QL, DESeq2, or limma-voom. There is no complementarity to exploit: on
this run the four methods order genes almost identically and differ only in
where they stop.

| Method | Calls | Extra over BARCS | Of which active | Precision of the extra calls |
|---|---:|---:|---:|---:|
| BARCS-original | 59 | --- | --- | --- |
| edgeR-QL | 93 | 34 | 16 | 0.47 |
| DESeq2 | 95 | 36 | 17 | 0.47 |
| limma-voom | 92 | 33 | 16 | 0.48 |

So "they have more false positives" is true but uninformative on its own: the
extra calls are roughly a coin flip, not noise. The informative version is
that the coin is not fair everywhere. Stratifying edgeR-QL's 34 extra calls by
how close BARCS came to calling them:

| BARCS gene FDR | Extra edgeR-QL calls | Active | Precision |
|---|---:|---:|---:|
| 0.10--0.15 | 8 | 5 | 0.63 |
| 0.15--0.25 | 14 | 9 | 0.64 |
| 0.25--0.50 | 11 | 2 | 0.18 |
| above 0.50 | 1 | 0 | 0.00 |

The extra discoveries separate cleanly into a recoverable band and a junk
tail. Genes sitting at BARCS FDR 0.10--0.25 are about 63% likely to be truly
active; past 0.25 the count models are mostly picking up noise. A BARCS screen
reporting genes in the 0.10--0.25 band as a clearly labelled secondary tier
would recover most of the difference without adopting the count models'
calibration.

**The 16 genes all three catch and BARCS misses are real but weak.** Their
median absolute simulated phenotype is 0.389, against 0.625 for the 55 actives
BARCS does call (0.581 across all 79 actives). Fourteen of the 16 sit at BARCS
FDR between 0.11 and 0.25 --- just over the line, not far from it. They are
not a distinct biological class: 12 linear and 4 sigmoidal, against 40 and 15
among the genes BARCS catches.

**Most of the gap is the calibration ratchet, not the count model.** On this
run `bb_calibrate_controls()` estimates a scale of 1.238, which is real work:
the uncalibrated negative-control $p<0.05$ rate is 0.08 rather than 0.05.
But sweeping the scale shows it overshoots at gene level:

| Control scale | Calls | True positives | Realized FDP | Missed genes recovered |
|---:|---:|---:|---:|---:|
| 1.00 (none) | 87 | 71 | 0.184 | 14 of 16 |
| 1.08 | 71 | 61 | 0.141 | 6 |
| 1.16 | 67 | 60 | **0.104** | 5 |
| 1.24 (applied) | 59 | 55 | 0.068 | 0 |
| 1.40 | 53 | 52 | 0.019 | 0 |

At the applied scale BARCS spends its whole FDR budget and then some: it lands
at realized FDP 0.068 when 0.10 was permitted. A scale of 1.16 would sit on the
nominal target and return five more true positives. Removing calibration
entirely recovers 14 of the 16 but pushes FDP to 0.184, so the ratchet is
buying something --- it is simply set tighter than the nominal threshold
requires. This is the clearest actionable finding in the comparison: the
`alpha = 0.05` control-calibration target and the 0.10 gene-FDR target are not
consistent with each other.

**Two genes resist even that.** GENE0063 and GENE0314 stay uncalled with no
calibration at all, and both fail for the guide-level reason:

- GENE0314 (phenotype $+0.672$) has two guides with essentially no knockdown
  (0.002 and 0.009) that contribute noise, and among its three real guides
  `g01570` carries a strong correct-signed effect ($\widehat\beta=+0.295$) that
  BARCS scores at $p=0.230$ because its own dispersion estimate came out high
  (variance inflation 10.3 against a library median of 2.7). edgeR-QL, which
  borrows dispersion across the library, scores the same guide at $p=0.080$.
- GENE0063 loses one guide to the `min_total_count = 30` filter and has one
  genuinely discordant guide, leaving three concordant guides against an
  equal-weight signed-$z$ combination.

These are the per-guide-dispersion mechanism made concrete: with seven
residual degrees of freedom the unshared estimate is noisy, and when it draws
high on an informative guide that guide's evidence is thrown away. Note that
this is a stochastic per-guide failure, not a systematic depth effect ---
across the library, variance inflation is uncorrelated with guide read depth
(Spearman $-0.03$, flat across read quartiles).

### Closing the gap: guide-dispersion moderation

The deficits above have one cause. BARCS estimates a beta-binomial dispersion
per guide from that guide's own residuals, on the $2R-1$ residual degrees of
freedom the design supplies. That estimate is noisy in both directions, and the
damaging direction is the quiet one: a few guides whose residuals happen to
look unusually clean report standard errors their data do not support. Those
guides drive the negative-control tail, and `bb_calibrate_controls()` can only
answer that tail with a single scale applied to *every* guide in the screen.

`bb_moderate_dispersion()` corrects those guides individually instead. It
shrinks each guide's variance inflation --- the untruncated Pearson dispersion
`pearson_null / df_residual` --- toward a trend fitted across the library, using
the scaled-F moment estimator of Smyth (2004). Both the prior degrees of
freedom and the shrinkage target are estimated from the screen, so there is no
tuning constant that could be fitted to an outcome.

The direction of the effect is worth stating plainly, because it is not the
obvious one. Moderation makes genuinely noisy guides *stricter*. What buys the
power is that the quiet guides stop inflating the empirical null, so the
blanket penalty drops: on the baseline realization the control scale falls from
1.238 to 1.172, and to 1.105 when the original tail-quantile estimator is kept.
More discoveries at a lower realized FDP is the net result of testing most
guides more sharply and a few more strictly.

A second, smaller change targets the null-scale estimator itself.
`bb_calibrate_controls()` matched one order statistic --- the empirical
$1-\alpha$ quantile of the absolute control statistics --- to the corresponding
$t$ cutoff. With a few hundred controls that is high-variance.
`method = "qq_slope"` fits the slope of the control quantile-quantile plot
against the $t$ reference over the 0.50--0.95 band instead, using many order
statistics rather than one.

#### How this was validated

Tuning a method on the benchmark it is reported on is how methods get oversold,
so the work was split three ways, with simulations regenerated from the pinned
Julia environment in `julia/` (re-simulating the committed baseline reproduces
`counts.tsv` byte for byte):

- **Development** --- the baseline scenario at the five seeds already used in
  this repository. Every design decision was made here.
- **Held-out** --- the nine supported multi-replicate scenarios at three seeds
  never inspected during development.
- **Confirmatory** --- the same nine scenarios at three further fresh seeds,
  generated and run only after the held-out split had overturned one default.

Held-out means over 27 runs, and confirmatory means over 26:

| Method | AP (held-out) | F1 (held-out) | FDP (held-out) | AP (confirm.) | F1 (confirm.) | FDP (confirm.) |
|---|---:|---:|---:|---:|---:|---:|
| BARCS-moderated | 0.8722 | **0.7852** | **0.0790** | **0.9077** | 0.8160 | **0.0852** |
| edgeR-QL | 0.8718 | 0.7837 | 0.2137 | 0.9051 | 0.8141 | 0.1990 |
| DESeq2 | **0.8724** | 0.7830 | 0.2315 | 0.9062 | 0.8113 | 0.2243 |
| limma-voom | 0.8699 | 0.7813 | 0.2177 | 0.9040 | 0.8094 | 0.2137 |
| BARCS-original | 0.8547 | 0.7289 | 0.0777 | 0.8898 | 0.7801 | 0.0803 |

Against BARCS-original the gain is unambiguous and repeats on both splits:
average precision $+0.018$ (95% interval $+0.013$ to $+0.022$ held-out,
$+0.014$ to $+0.022$ confirmatory) and F1 $+0.056$ and $+0.036$, at no
measurable cost in realized FDP.

Against the count models the honest statement is a **tie on ranking and F1,
won on calibration**. Paired F1 differences are within noise on both splits
(largest interval $\pm0.025$), average-precision differences are at most
$0.004$ and change sign between splits, and the realized-FDP difference is
$0.11$ to $0.15$ in BARCS's favour on every comparison, always significant.
BARCS now matches the count models where it used to trail them, and keeps a
false-discovery proportion under the nominal 0.10 target while theirs sits at
two to three times nominal. That is the result the earlier sections were
missing.

#### Two things the splits caught

Neither would have been visible from a single evaluation, and both are recorded
rather than folded away.

- **A default chosen on development data was wrong.** The development split
  favoured a conservative moderation --- one-way, so no guide variance is ever
  lowered, and without claiming the prior degrees of freedom --- because
  textbook two-sided moderation ran at FDP 0.163 on the baseline scenario. On
  both fresh splits that caution proved unnecessary and expensive: two-sided
  moderation held FDP at 0.079 and 0.085, while the conservative variant gave
  up F1 (paired $-0.033$ held-out, $-0.020$ confirmatory, both significant).
  The default is therefore the standard two-sided form, with `one_way` and
  `borrow_df` retained as an explicitly conservative option. The confirmatory
  split exists because that reversal was decided on held-out data.
- **An empirical-null idea was discarded, not shipped.** Rescaling the gene
  statistic against a central-window normal fit (Efron's empirical null) looked
  strong on the baseline realization. It is not in the codebase: on synthetic
  checks it is biased even with no active genes, and with one-sided activity
  --- ordinary in a dropout screen --- it returned a scale of 10.97 against a
  truth of 1.5. Its apparent success depended on this simulator's roughly
  symmetric mix of increasing and decreasing genes.

One caveat: three held-out scenarios (MOI 0.10, MOI 0.40, and the 100-gene
library) run above nominal FDP for BARCS both before and after moderation, so
the calibration claim is a screen-average one rather than a per-scenario
guarantee. One confirmatory run is skipped entirely because a 100-gene library
left fewer than the 20 negative-control guides `bb_calibrate_controls()`
requires; that limit predates the moderation and stops BARCS-original equally.

### Genome scale: is FDR control robust across cutoffs?

A single realized-FDP number cannot distinguish a method that controls the
false-discovery rate from one that happens to land near 0.10 — and 0.10 is
looser than a genome-wide screen would use anyway. So the check is repeated on
a log-spaced grid of thirteen cutoffs from 1e-6 to 0.20
(`examples/crispulator_facs_moi_10k_benchmark.R`, MOI 0.20, three seeds,
common gene set). The quantity is realized FDP ÷ requested cutoff: one is
exact, above one means more false discoveries than advertised.

| Method | 1e-6 | 1e-5 | 1e-4 | 1e-3 | 0.01 | 0.05 | 0.10 | Cells held |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BARCS-original | 0.0 | 0.0 | 0.0 | 0.6 | 0.6 | 0.6 | 0.7 | 37/39 |
| BARCS-moderated | 0.0 | 0.0 | 0.0 | 0.4 | 0.7 | 0.6 | 0.7 | 37/39 |
| MAGeCK-MLE | 0.0 | 0.0 | 0.0 | 0.0 | 0.1 | 0.1 | 0.1 | **39/39** |
| CRISPhieRmix | 0.0 | 0.0 | 0.0 | 0.0 | 1.7 | 2.0 | 2.0 | 22/39 |

"Cells held" = cutoff-by-seed combinations, of 13 × 3 = 39, where realized FDP
did not exceed the cutoff.

**CRISPhieRmix overshoots about twofold wherever it calls enough to measure**
(22/39 cells). Its zeros at the strict end are not calibration — it returns 0
and 4 genes at 1e-6 and 1e-5, so there is nothing there to be wrong about.
That strict end is precisely the regime a genome-wide screen is read in.

Both BARCS statistics track the requested rate down to 1e-4 and return no
false discoveries below it (37/39 cells), while still calling 276 genes at
1e-6. MAGeCK-MLE holds all 39 but never spends more than a tenth of its
budget.

The general count models behaved the same way but far more extremely; that
analysis is retained in the 400-gene section above, and it is what first
established that a fixed nominal cutoff means different things to different
methods.

**BARCS exceptions, named not rounded away.** Both miss 2 of 39 cells, all
from one seed (20250725) at the 3e-4 and 1e-3 cutoffs, each a single false
discovery against an allowance below one gene. Monte Carlo granularity, not
systematic failure — and the reason the calibration claim is a screen average
rather than a finite-sample guarantee. At MOI 0.30, BARCS-moderated holds
35/39.

#### F1 across cutoffs

| Method | 1e-6 | 1e-5 | 1e-4 | 1e-3 | 0.01 | 0.05 | 0.10 |
|---|---:|---:|---:|---:|---:|---:|---:|
| BARCS-original | 0.046 | 0.126 | 0.262 | 0.456 | 0.653 | 0.774 | 0.811 |
| BARCS-moderated | 0.248 | 0.354 | 0.479 | 0.607 | 0.748 | 0.827 | **0.847** |
| MAGeCK-MLE | **0.306** | 0.306 | 0.306 | 0.306 | 0.403 | 0.622 | 0.718 |
| CRISPhieRmix | 0.000 | 0.004 | 0.080 | 0.321 | 0.634 | 0.780 | 0.785 |

Restricted to CRISPR-specific callers the ordering is stable: BARCS-moderated
leads at every cutoff from 1e-5 to 0.10. The lone exception is 1e-6, where
MAGeCK-MLE's 0.306 is a quantization artifact — its call count is identically
358 from 1e-6 to 3e-3 — rather than a measurement. CRISPhieRmix is the mirror
image, collapsing to 0.000 and 0.004 at the strict end because it has stopped
calling.

That stability is a consequence of the comparison set, and worth stating as
such: a fixed requested cutoff does not place different methods at the same
real error rate, so a method that overshoots looks strong on F1 for that
reason alone.

#### Matched on realized FDP

| Method | 0.001 | 0.002 | 0.005 | 0.01 | 0.02 | 0.05 | 0.10 |
|---|---:|---:|---:|---:|---:|---:|---:|
| BARCS-original | 0.491 | 0.529 | 0.631 | 0.678 | 0.738 | 0.795 | 0.819 |
| BARCS-moderated | 0.640 | 0.677 | 0.724 | 0.764 | 0.808 | 0.838 | 0.847 |
| MAGeCK-MLE | 0.469 | 0.509 | 0.666 | 0.739 | 0.789 | — | — |
| CRISPhieRmix | 0.372 | 0.415 | 0.488 | 0.553 | 0.640 | 0.725 | 0.775 |
| **Spread** | 0.268 | 0.262 | 0.236 | 0.211 | 0.168 | 0.113 | 0.072 |

**BARCS-moderated leads at every matched realized FDP**, from 0.640 against
0.372–0.491 at a matched 0.001 to 0.847 against 0.775–0.819 at 0.10. Holding
the error rate fixed is the point: at the same realized error rate, moderated
BARCS recovers more genes than either CRISPR-specific comparator. The lead
narrows as the budget grows, 0.268 of spread at 0.001 down to 0.072 at 0.10,
which is what should happen when every method is allowed to be wrong more
often.

Two caveats. MAGeCK-MLE's two strictest entries rest on its quantized FDR and
are indicative only, and it has no entry at matched 0.05 or 0.10 because it
never becomes that liberal — absence there is conservatism, not failure.
**BARCS-MOD over BARCS-ST is the cleanest comparison**, larger at every
matched error rate, so moderation is a real gain within the BARCS family
rather than a threshold artifact.

Note that this reverses the earlier reading. When the general count models
were in the comparison, F1 converged above a matched 0.005 and no method
separated. Restricted to CRISPR-specific callers the separation is clear
throughout — the count models were the methods that closed the gap, by
operating at error rates BARCS never reaches.

#### MAGeCK-MLE is not resolvable at strict cutoffs

MAGeCK-MLE's call count is constant across orders of magnitude of requested
FDR. That is a property of its inference, not of the data: its gene p-value is
a permutation tail probability from `genes × rounds` draws, quantized in steps
of ~1/(genes × rounds), and any gene whose statistic beats every permuted
value gets p exactly 0 — so it is called at *any* cutoff.

Refitting one realization at 10 rounds instead of 1
(`examples/crispulator_facs_mageck_permutation_resolution.R`) confirms the
mechanism and shows it is **not** a settings artifact that can be corrected:

| Rounds | Distinct p | Smallest positive p | Predicted 1/(genes×rounds) | Genes at p = 0 |
|---:|---:|---:|---:|---:|
| 1 | 3,527 | 2e-4 | 1.06e-4 | 186 |
| 10 | 7,713 | 2e-5 | 1.06e-5 | 66 |

The resolution improves exactly as predicted. But the plateau persists — and
moves *down*:

| Requested FDR | Calls (1 rd) | Calls (10 rd) | F1 (1 rd) | F1 (10 rd) |
|---|---:|---:|---:|---:|
| 1e-6 … 1e-3 | 186 | 66 | 0.174 | 0.065 |
| 0.01 | 453 | 488 | 0.376 | 0.399 |
| 0.10 | 1,171 | 1,141 | 0.744 | 0.732 |

More permutation lowers the plateau rather than removing it, because fewer
genes reach p = 0. MAGeCK-MLE's strict-cutoff position is therefore set by a
computational budget, not by the screen. Resolving FDR 1e-6 at this library
size needs ~1e6/9475 ≈ 110 rounds, about fourteen hours per realization at the
throughput measured here.

**We therefore read MAGeCK-MLE's curve below a requested 0.01 as not
resolvable rather than as a performance result**, and exclude it from the
strict-regime comparison. The same caution applies to any permutation-based
gene FDR at genome scale.

### The one-replicate boundary, and why the count models lead there

At $R=1$ the design has three samples (low, bulk, high) and the guide model
drops to `~ phenotype_z`, leaving **one residual degree of freedom**. This
setting is excluded from the headline comparison above; the numbers below are
the `one-replicate boundary` scope, five seeds.

| Method | Average precision | AUROC | Directional recall | Realized FDP | F1 |
|---|---:|---:|---:|---:|---:|
| BARCS-original | 0.633 | 0.790 | 0.000 | 0.000 | 0.000 |
| BARCS-EB | 0.632 | 0.792 | 0.117 | 0.000 | 0.201 |
| MAGeCK-MLE | 0.548 | 0.780 | 0.014 | 0.000 | 0.026 |
| edgeR-QL | 0.696 | 0.815 | 0.708 | 0.526 | 0.564 |
| DESeq2 | 0.694 | 0.814 | 0.419 | 0.141 | 0.554 |
| limma-voom | 0.698 | 0.816 | 0.684 | 0.507 | 0.571 |

Two different things are happening, and only one of them is a genuine
advantage.

**The F1 and recall lead is mostly not real.** edgeR-QL calls with realized
FDP 0.526 and limma-voom with 0.507: more than half of their discoveries are
false. Their negative-control $p<0.05$ rates are 0.228 and 0.209, four times
nominal. An F1 of 0.564 purchased at FDP 0.526 is not a usable operating
point, it is a method that has stopped controlling anything. DESeq2 is the
honest exception, holding FDP to 0.141, and it is the only count model whose
$R=1$ calls are defensible.

**The ranking lead is real, and it is instructive.** Average precision is
0.696 for edgeR-QL against 0.633 for BARCS-original, a paired difference of
$+0.063$ (95% interval $+0.030$ to $+0.096$) --- nearly four times the
$+0.017$ gap in the multi-replicate settings. Two mechanisms compound:

1. *The dispersion estimate collapses.* BARCS estimates each guide's
   overdispersion from that guide's own residuals. The relative variance of
   such an estimate scales roughly as $2/d$ in the residual degrees of
   freedom $d$: about 0.29 at the four-replicate baseline ($d=7$), but about
   2.0 at $R=1$ ($d=1$). The standard error is then nearly pure noise, and
   the ranking it induces degrades with it. edgeR-QL, limma-voom, and DESeq2
   borrow dispersion from the whole library, so their per-guide variance is
   dominated by the prior and barely moves as $d$ falls. Information sharing
   is merely helpful at $d=7$; at $d=1$ it is what keeps the ranking alive.
2. *The $t_1$ reference compresses the evidence.* With one residual degree of
   freedom the guide reference distribution is Cauchy, so even $|t|=10$ gives
   $p\approx0.06$. No guide can contribute a small $p$-value, the dynamic
   range of the guide evidence is squashed into a narrow band before the
   signed-$z$ combination, and relative differences between strong and weak
   guides are flattened. The transform is monotone within a guide but the
   Stouffer sum across guides is not, so gene-level ranking is genuinely
   lost. This is the effect that motivates the separate BARCS-GC
   guide-consistency statistic.

The zero calls follow from the same place: a one-degree-of-freedom $t$ test
cannot produce $p$-values small enough to survive Benjamini--Hochberg across
400 genes. BARCS-GC is the prespecified response, raising average precision to
0.670, directional recall to 0.145, and F1 to 0.245 at realized FDP 0. Its
calls are hypothesis-generating rather than confirmatory.

The practical reading is that $R=1$ is not a regime where the count models are
better calibrated --- two of the three are far worse. It is a regime where
borrowing variance across guides is the only way to retain a usable gene
ranking, and where no method in this comparison should be trusted to control
FDR.

### What is being compared at gene level?

MAGeCK-MLE supplies its native joint gene-level beta and Wald inference.
edgeR-QL, DESeq2, and limma-voom are fitted guide by guide and then use the
same signed-$z$ guide aggregation as BARCS-original. This makes their
guide-level count models comparable while holding the historical gene
combiner fixed, but it is not a claim that signed-\(z\) aggregation is the
only or canonical gene-level interface for those packages.

The four BARCS methods reuse one shared beta-binomial guide fit:
BARCS-original combines calibrated signed guide scores; BARCS-NORM estimates
the arithmetic mean and sample standard deviation of guide beta values;
BARCS-partial fits a measurement-error random-effects gene mean; and BARCS-EB
moderates guide heterogeneity toward an empirical prior.

## Waterbear GSE242880

This real four-bin IL2RA screen has 26 directionally validated regulators and
seven additional tested candidates that did not validate. BARCS and MAGeCK
were evaluated from complete per-gene outputs. Waterbear and MAUDE are
literature-reported aggregates only.

| Method | Design | Screen calls | Directionally validated | Selected-panel F1 | Average precision |
|---|---|---:|---:|---:|---:|
| BARCS-original | four-bin trend + donor | 49 | 22/26 | 0.863 | **0.945** |
| BARCS-NORM | four-bin trend + donor | 0 | 0/26 | 0.000 | 0.897 |
| BARCS-partial | four-bin trend + donor | 71 | 19/26 | 0.792 | 0.910 |
| BARCS-EB | four-bin trend + donor | 60 | **23/26** | **0.885** | 0.940 |
| MAGeCK-MLE | four-bin trend + donor | 72 | 17/26 | 0.723 | 0.872 |
| MAGeCK test | outer Q1 versus Q4 | 30 | 18/26 | 0.766 | 0.906 |
| Waterbear | published four-bin model | 79 | 24/26 | not available | not available |
| MAUDE | published four-bin model | 406 | 25/26 | not available | not available |

BARCS-EB is the best of the six rerun analyses for validated recovery and
selected-panel F1. BARCS-original makes fewer calls and has the best selected-
panel ranking. BARCS-NORM retains useful ranking but makes no FDR-0.10 calls,
which is consistent with its small-guide t-reference penalty. The reported
Waterbear and MAUDE recoveries are higher, but
they cannot be assigned F1 or average precision without their complete
per-gene results for the seven non-validating candidates. MAUDE's 25/26
recovery also accompanies 406 calls, so it is sensitivity rather than a
complete specificity result.

This dataset favors specialist joint-bin modeling on likelihood fidelity:
four FACS bins from one donor partition one pool and are not independent.
BARCS remains useful as a transparent ordered-trend analysis with explicit
negative-control diagnostics, but it does not replace Waterbear's joint
multinomial model.

## Liang Cas13 processed-count sensitivity

The Liang comparison uses rounded versions of the deposited normalized,
ComBat-corrected Day-0 and Day-14 values. All newly fitted methods receive
that same matrix and replicate/day design. This controls the input comparison,
but it is not a raw-count likelihood benchmark.

At gene level, BARCS-original, edgeR-QL, DESeq2, and limma-voom use the same
unweighted signed-\(z\) rule:

$$
Z_g=\frac{1}{\sqrt{m_g}}\sum_{j=1}^{m_g}
\operatorname{sign}(\widehat\beta_{gj})
\Phi^{-1}\!\left(1-\frac{p_{gj}}{2}\right).
$$

Thus these four methods differ in their guide-level count model, not in the
final gene combiner. MAGeCK-MLE uses its native joint gene model.
BARCS-partial and BARCS-EB are different again: they pool guide effects by
random-effects inverse-variance weighting and do not use Stouffer
aggregation.

In HAP1, the extra guide-level models recover more of the 60 essential
controls at FDR 0.10, while BARCS has the best nominal null specificity:

| Method | Essential controls recovered | Null specificity at \(p<0.05\) |
|---|---:|---:|
| BARCS | 44/60 | **0.957** |
| MAGeCK-MLE | 49/60 | 0.895 |
| edgeR-QL | **52/60** | 0.931 |
| DESeq2 | **52/60** | 0.941 |
| limma-voom | **52/60** | 0.929 |

The primary analysis retains every valid guide. A pre-outcome sensitivity
selecting the five guides with greatest mean Day-0 abundance reduces BARCS
macro-average precision from 0.776 to 0.638 and FDR-0.10 essential recall
from 0.487 to 0.330. Selecting guides by their observed p-values would be
circular and is not used.

## Reproducibility and provenance

The CRISPulator metrics were recomputed from all method-level gene-result
files with one evaluator. The external fits were present from a previous
complete run; they were not silently described as newly refitted. Before
using them, the comparison verifies that the historical BARCS effect vector
matches the current BARCS-original vector exactly in every run and that
$p$-value/FDR differences are below $10^{-3}$. The largest observed
$p$-value difference is $1.95\times10^{-4}$, caused by a small
negative-control scale change; effect differences are zero.

The provenance table records the MD5 hashes of every count matrix, design,
truth table, and result file. Thus the compact comparison can be audited even
when large ignored result directories are not committed.

Each analysis is split into a fitting stage, which needs the ignored
`results/` tree, and an aggregation stage, which reads only the committed
per-run table. The aggregation stage owns the scope rule and can be re-run on
its own to rebuild the summary tables and figures:

```
Rscript examples/crispulator_facs_external_head_to_head_aggregate.R
Rscript examples/crispulator_facs_f1_threshold_curves_aggregate.R
```

- `examples/crispulator_facs_external_head_to_head.R`
- `examples/crispulator_facs_external_head_to_head_aggregate.R`
- `data/derived/crispulator_facs_external_head_to_head_metrics.csv`
- `data/derived/crispulator_facs_external_head_to_head_provenance.csv`
- `examples/waterbear_facs_external_head_to_head.R`
- `data/derived/waterbear_facs_external_head_to_head_metrics.csv`
- `examples/crispulator_facs_f1_threshold_curves.R`
- `examples/crispulator_facs_f1_threshold_curves_aggregate.R`
- `data/derived/crispulator_facs_f1_by_fdr.csv`
- `figures/crispulator_facs_f1_by_fdr.pdf`

The moderation and its validation:

- `examples/crispulator_facs_improved_barcs.R` (committed inputs only)
- `data/derived/crispulator_facs_improved_barcs_baseline.csv`
- `data/derived/crispulator_facs_improved_barcs_scales.csv`
- `examples/crispulator_facs_improved_barcs_holdout.R` (`--simulate` first;
  needs Julia and the pinned environment in `julia/`)
- `data/derived/crispulator_facs_improved_barcs_holdout_metrics.csv`
- `data/derived/crispulator_facs_improved_barcs_holdout_summary.csv`
- `data/derived/crispulator_facs_improved_barcs_holdout_paired.csv`

The per-gene case study runs entirely from committed inputs and needs
edgeR, limma, and DESeq2:

- `examples/crispulator_facs_miss_cases.R`
- `data/derived/crispulator_facs_miss_cases.csv`
- `data/derived/crispulator_facs_miss_marginal_value.csv`
- `data/derived/crispulator_facs_miss_precision_bands.csv`
- `data/derived/crispulator_facs_calibration_scale_sweep.csv`
- `examples/liang_hap1_specificity_volcano.R`
- `figures/liang_hap1_specificity_volcano.pdf`
