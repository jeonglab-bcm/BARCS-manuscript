# BARCS compared with other methods

This report separates two questions that can otherwise give contradictory
answers:

1. **Ranking:** are truly active genes near the top?
2. **Decision calibration:** when a method reports FDR below 0.10, is the
   realized false-discovery proportion actually near or below 0.10?

No method wins both questions in every benchmark.

## CRISPulator FACS simulations

The comparison uses the same 50 primary low + bulk + high simulations as the
three-method BARCS benchmark: ten one-at-a-time parameter scenarios and five
fixed seeds. Every method receives the same guide counts, ordered-phenotype
score, replicate adjustment, and gene truth.

| Method | Average precision | AUROC | Directional recall at FDR 0.10 | Realized FDP | F1 | Negative-control $p<0.05$ |
|---|---:|---:|---:|---:|---:|---:|
| BARCS-original | 0.856 | 0.918 | 0.638 | 0.085 | 0.702 | **0.055** |
| BARCS-partial | 0.779 | 0.885 | 0.390 | 0.052 | 0.525 | 0.016 |
| BARCS-EB | 0.840 | 0.910 | 0.463 | **0.013** | 0.599 | 0.013 |
| MAGeCK-MLE | 0.849 | 0.922 | 0.587 | 0.058 | 0.683 | 0.037 |
| edgeR-QL | **0.877** | **0.928** | **0.822** | 0.277 | 0.765 | 0.138 |
| DESeq2 | 0.877 | 0.927 | 0.795 | 0.233 | 0.769 | 0.112 |
| limma-voom | 0.877 | 0.928 | 0.813 | 0.257 | **0.772** | 0.128 |

Bold identifies the largest ranking, recall, or F1 value, the smallest
realized FDP, and the negative-control rate closest to its nominal 0.05
target. A very small FDP is conservative rather than automatically optimal;
it should be read together with recall.

The result is a real trade-off:

- edgeR-QL, DESeq2, and limma-voom rank simulated active genes slightly
  better and recover many more at their nominal threshold.
- Their nominal FDR 0.10 calls are anti-conservative here: mean realized FDP
  is 0.23--0.28, and 0.11--0.14 of negative controls have \(p<0.05\).
- BARCS-original is the strongest BARCS ranker and its nominal threshold is
  much closer to the intended operating point. Its average precision is
  about 0.022 below edgeR-QL, but its realized FDP is about 0.192 lower.
- BARCS-EB is the safest method, but the cost is substantial recall.
- MAGeCK-MLE is well calibrated and close to BARCS-original in ranking. On
  paired runs, its average-precision difference from BARCS-original is
  $-0.0065$, with a 95% interval from $-0.0163$ to $0.0033$; its F1 is
  lower by 0.0191 (95% interval $-0.0345$ to $-0.0036$).

### F1 across nominal FDR thresholds

For a called gene set, the standard binary F1 score is

$$
\mathrm{F1}=\frac{2TP}{2TP+FP+FN}.
$$

Here a true positive is any simulated active gene passing the nominated
gene-FDR threshold; direction is reported separately. The threshold scan uses
0.10, 0.05, 0.01, 0.005, and 0.001 in every one of the 50 runs.

The scan changes the interpretation. At nominal FDR 0.10, limma-voom has the
largest mean F1 (0.772) but realized FDP 0.257. At nominal FDR 0.01, edgeR-QL
has mean F1 0.796 and realized FDP 0.087. BARCS-original at nominal FDR 0.10
has F1 0.702 and realized FDP 0.085. Thus the general count models retain a
real power advantage after moving to a threshold with comparable empirical
FDP in this simulator; their advantage is not explained entirely by the
anti-conservative 0.10 operating point.

This is a diagnostic, not permission to select a method-specific threshold
from the known truth in new data. The selected 0.01 threshold would need
independent calibration or negative controls before prospective use.

At the four-replicate baseline, edgeR-QL has the highest mean average
precision (0.932), whereas BARCS-original has 0.917. limma-voom has the
highest F1 (0.832) and BARCS-original has 0.828, but their realized FDPs are
0.218 and 0.086, respectively.

### What is being compared at gene level?

MAGeCK-MLE supplies its native joint gene-level beta and Wald inference.
edgeR-QL, DESeq2, and limma-voom are fitted guide by guide and then use the
same signed-$z$ guide aggregation as BARCS-original. This makes their
guide-level count models comparable while holding the historical gene
combiner fixed, but it is not a claim that signed-\(z\) aggregation is the
only or canonical gene-level interface for those packages.

The three BARCS methods reuse one shared beta-binomial guide fit:
BARCS-original combines calibrated signed guide scores, BARCS-partial fits a
random-effects gene mean, and BARCS-EB moderates guide heterogeneity toward an
empirical prior.

## Waterbear GSE242880

This real four-bin IL2RA screen has 26 directionally validated regulators and
seven additional tested candidates that did not validate. BARCS and MAGeCK
were evaluated from complete per-gene outputs. Waterbear and MAUDE are
literature-reported aggregates only.

| Method | Design | Screen calls | Directionally validated | Selected-panel F1 | Average precision |
|---|---|---:|---:|---:|---:|
| BARCS-original | four-bin trend + donor | 49 | 22/26 | 0.863 | **0.945** |
| BARCS-partial | four-bin trend + donor | 71 | 19/26 | 0.792 | 0.910 |
| BARCS-EB | four-bin trend + donor | 60 | **23/26** | **0.885** | 0.940 |
| MAGeCK-MLE | four-bin trend + donor | 72 | 17/26 | 0.723 | 0.872 |
| MAGeCK test | outer Q1 versus Q4 | 30 | 18/26 | 0.766 | 0.906 |
| Waterbear | published four-bin model | 79 | 24/26 | not available | not available |
| MAUDE | published four-bin model | 406 | 25/26 | not available | not available |

BARCS-EB is the best of the five rerun analyses for validated recovery and
selected-panel F1. BARCS-original makes fewer calls and has the best selected-
panel ranking. The reported Waterbear and MAUDE recoveries are higher, but
they cannot be assigned F1 or average precision without their complete
per-gene results for the seven non-validating candidates. MAUDE's 25/26
recovery also accompanies 406 calls, so it is sensitivity rather than a
complete specificity result.

This dataset favors specialist joint-bin modeling on likelihood fidelity:
four FACS bins from one donor partition one pool and are not independent.
BARCS remains useful as a transparent ordered-trend analysis with explicit
negative-control diagnostics, but it does not replace Waterbear's joint
multinomial model.

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

- `examples/crispulator_facs_external_head_to_head.R`
- `data/derived/crispulator_facs_external_head_to_head_metrics.csv`
- `data/derived/crispulator_facs_external_head_to_head_provenance.csv`
- `examples/waterbear_facs_external_head_to_head.R`
- `data/derived/waterbear_facs_external_head_to_head_metrics.csv`
- `examples/crispulator_facs_f1_threshold_curves.R`
- `data/derived/crispulator_facs_f1_by_fdr.csv`
- `figures/crispulator_facs_f1_by_fdr.pdf`
- `examples/liang_hap1_specificity_volcano.R`
- `figures/liang_hap1_specificity_volcano.pdf`
