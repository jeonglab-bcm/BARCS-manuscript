# BARCS re-analysis of the Liu et al. in vivo T cell CRISPR screen

**Prepared by:** Kelly Lee (KyuWon.Lee@bcm.edu), Liuzlab-DSC
**Date:** 24 August 2026
**Source:** Liu Q, Chen PA, Urs E, *et al.* "In vivo genome-wide CRISPR screens of human T cells in solid tumours." *Nature* (2026). doi:10.1038/s41586-026-10906-9 — Supplementary Table 5.

---

## One-paragraph summary

I refit the focused sub-library screen from Liu et al. using BARCS, a guide-level
beta-binomial regression that keeps donor and mouse in the design rather than
collapsing each contrast to treatment-vs-control. The fit is well-behaved — the 52
non-targeting controls come out null, the empirical-null calibration factor is 1.00
(no correction needed), 100% of guides converged, and the screen's built-in positive
controls IFNG and TNF are the two strongest signals in the correct direction. BARCS
recovers the published hits, including GNAS. The one result the original analysis
could not produce comes from pairing on the animal instead of the donor: **P2RY8 and
PTGER4 do not survive that pairing**, indicating their apparent signal was
between-mouse variance rather than biology. This is a re-analysis of published
normalised counts, not an independent replication.

## What was done

| Item | Description |
|---|---|
| Data | Supplementary Table 5 — the only screen in the paper publishing per-mouse counts |
| Library | 268 sgRNAs · 36 genes × 6 guides · 52 non-targeting controls (NTCTRL) |
| Readout | Tumour-infiltrating T cells sorted by IFNγ production, top 20% vs bottom 20% |
| Arm A | CD3-scFv platform vs A375low — 16 libraries, 8 mice, 3 donors |
| Arm B | NY-ESO-1 TCR vs WT A375 — 8 libraries, 4 mice, 2 donors |
| Models | `~ gate + donor` on each arm; `~ gate + mouse` on Arm A as a sensitivity check |

A positive coefficient means the knockout is enriched among IFNγ-high T cells — i.e.
losing that gene makes T cells better effectors.

**This is not the model the paper used.** The published analysis called hits with
MAGeCK, which collapses each contrast to a single treatment-vs-control comparison and
does not carry donor or mouse as explicit terms. BARCS instead fits a guide-level
beta-binomial regression with those factors kept in the design, so the two are
comparable in the **sign and rank** of gene effects but not in their numeric
magnitudes (see Limits #5). Keeping donor and mouse in the model is exactly what makes
the Result 3 pairing analysis possible — it is the intended point of difference, not
an attempt to reproduce MAGeCK's arithmetic.

## Result 1 — the fit is trustworthy

| Check | Arm A (`~ gate + donor`) |
|---|---|
| Guides with finite fits | 267 / 268 |
| IRLS convergence | 100% |
| Dispersion-boundary hits | 0% |
| Median extra-binomial ρ | 1.22 × 10⁻⁴ |
| Negative-control calibration factor | **× 1.00** |
| NTCTRL as a gene | effect −0.019, p = 0.26, FDR 0.63 — null |
| IFNG (positive control) | −0.538, FDR 1.6 × 10⁻⁴, 6/6 guides agree |
| TNF (positive control) | −0.266, FDR 2.1 × 10⁻⁴, 6/6 guides agree |

A calibration factor of exactly 1.00 means the non-targeting controls already behaved
as a proper null and BARCS found nothing to correct — the strongest single indicator
that the fit is sound.

## Result 2 — the published biology reproduces

Arm A hits at FDR < 0.10, all positive (knockout raises IFNγ):
**MED12 0.008 · STT3B 0.008 · TNFAIP3 0.008 · NFKBIA 0.074 · GNAS 0.087 · RASA2 0.098**

Across the two independent in vivo arms, gene-level effects correlate at
**Pearson r = 0.72** (p = 6 × 10⁻⁷, 36 genes), with **9/9 sign agreement** among genes
reaching FDR < 0.10 in either arm. Both targets the authors prioritised reappear —
GNAS in Arm A (0.087), STUB1 in Arm B (0.048) — each in the arm where the paper's own
analysis found it strongest.

## Result 3 — the new finding

Both sort gates come from the same animal, so mouse is the natural pairing factor.
Refitting Arm A as `~ gate + mouse` costs degrees of freedom (12 → 7) but lowers median
ρ to 8.42 × 10⁻⁵ and raises guide-level discoveries from 7 to 21.

| Gene | Effect | `+ donor` | `+ mouse` | Arm B | Reading |
|---|---|---|---|---|---|
| TNFAIP3 | +0.228 | 0.008 | **0.001** | 0.004 | Strengthens |
| RASA2 | +0.262 | 0.098 | **0.025** | 8.1e−4 | Strengthens |
| NFKBIA | +0.251 | 0.074 | 0.031 | 0.038 | Holds in both arms |
| MED12 | +0.246 | 0.008 | 0.036 | 0.647 | Holds in Arm A |
| STT3B | +0.247 | 0.008 | 0.025 | 0.487 | Holds in Arm A |
| GNAS | +0.204 | 0.087 | 0.105 | 0.208 | Borderline both ways |
| **P2RY8** | −0.105 | 0.109 | **0.555** | 0.881 | **Collapses** |
| **PTGER4** | +0.149 | 0.109 | **0.496** | 0.881 | **Collapses** |
| NTCTRL | −0.019 | 0.631 | 0.652 | 0.383 | Null, as required |

**Recommendation:** treat P2RY8 and PTGER4 as unsupported by this screen. Their signal
under the donor model is consistent with animal-to-animal bottleneck variance.

## Limits — state these before any number is quoted

1. **Re-analysis, not replication.** Every count came from the paper's own supplementary
   tables. Agreement with the published calls is expected and is not independent evidence.
2. **Counts are MAGeCK median-normalised and non-integer**, rounded to integers for the
   fit. The screen FASTQs were never deposited (only RNA-seq, GSE330227), so raw counts
   are not publicly obtainable. Library totals are uniform by construction and are not
   true sequencing depths; ranks and signs are sound, absolute dispersion is optimistic.
3. **The two arms were normalised in separate MAGeCK runs.** In Supplementary Table 5
   the CD3-scFv arm (Arm A) and the NY-ESO-1 TCR arm (Arm B) are each reported as
   already-normalised counts, and the median normalisation that produced them was run
   independently within each arm — so a guide's count in Arm A and the same guide's
   count in Arm B sit on different, non-comparable scales. Because the normalisation
   baselines differ, the arms cannot be placed in one model; they were fitted
   separately and their results compared only in sign and rank (this is why the r = 0.72
   cross-arm correlation in Result 2 is a concordance check, not a joint fit). A pooled
   `gate × system` interaction — the more powerful design — is not available without the
   raw per-library counts, which would let both arms be normalised together on one
   scale.
4. **GNAS is not settled by this table.** 0.087 with donor, 0.105 with mouse.
   Directionally consistent in all three fits, decisive in none. GNAS matters because it
   is one of the two targets the authors prioritised out of the whole screen: it encodes
   the stimulatory G-protein α-subunit (Gαs) that drives cAMP signalling, a pathway long
   argued to suppress T cell effector function, so a GNAS knockout that raises IFNγ is a
   mechanistically clean and therapeutically attractive hit — which is why the paper
   singles it out. This table gives GNAS the right sign in every fit but never crosses a
   confident FDR threshold, so it supports the paper's direction without independently
   confirming the hit. The paper's own case for GNAS rests on the genome-wide IFNγ screen
   plus individual-knockout validation, and nothing here contradicts that.
5. Effects are logit-scale beta-binomial coefficients — comparable to the published
   MAGeCK results in **sign and rank only**, not in magnitude.

## Reproducibility

Every output file is regenerated directly from the committed model inputs by
`examples/liu_tcell_barcs.R`, using the BARCS R package only (no BARCS Studio):

```
Rscript examples/liu_tcell_barcs.R
```

The script reads `results/liu_tcell/input/` and rewrites `results/liu_tcell/output/`.
For each arm it fits `bb_screen(term = "gate")`, applies non-targeting-control tail
calibration with `bb_calibrate_controls()`, and produces the gene-level table with
`bb_gene_consistency()`. BARCS is deterministic on this data; gene effects, the
calibration scale, and the hit list reproduce the numbers quoted above. One SPTLC2
guide sits on the beta-binomial convergence boundary and may be counted as converged
or not depending on the local linear-algebra backend, shifting a few guide-level FDRs
by ≈0.004 without changing any gene call.

## Suggested next step

Request the raw per-library sgRNA counts from the authors. That would allow (a) true
full-library denominators, (b) a pooled two-arm model with a `gate × system`
interaction, and (c) applying the same donor/mouse-aware design to the genome-wide
IFNγ screen, where the published analysis has the same unmodelled-variance issue and
far more genes at stake.

## Files

Generator: `examples/liu_tcell_barcs.R`.

Model inputs (`results/liu_tcell/input/`, regenerated from the published archive):
`armA_counts.csv`, `armA_metadata.csv`, `armB_counts.csv`, `armB_metadata.csv`

BARCS outputs (`results/liu_tcell/output/`; gene level 37 × 13 cols, guide level 268 × 17 cols):
`armA_gate-donor_GENES.csv` / `_GUIDES.csv`, `armA_gate-mouse_GENES.csv` / `_GUIDES.csv`,
`armB_gate-donor_GENES.csv` / `_GUIDES.csv`
