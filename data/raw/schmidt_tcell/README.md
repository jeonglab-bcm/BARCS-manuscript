# Schmidt et al. primary human T cell CRISPRa/CRISPRi screen sources

Source for the BARCS-versus-MAGeCK comparison in `results/schmidt_tcell/`.
Schmidt R, Steinhart Z, Layeghi M, *et al.*, "CRISPR activation and interference
screens decode stimulation responses in primary human T cells," *Science* 376,
eabj4008 (2022),
[doi:10.1126/science.abj4008](https://doi.org/10.1126/science.abj4008).

Nothing in this directory is versioned. One script fetches and prepares it:

```sh
Rscript scripts/prepare_schmidt_tcell.R
```

It needs the `readxl` package, downloads only what is missing, and writes the
eight model inputs into `results/schmidt_tcell/input/`.

## Downloaded series supplement

- `GSE174255_sgRNA-Read-Counts.xlsx` (22 MB) — the series supplementary file of
  GEO [GSE174255](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174255),
  `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174255/suppl/GSE174255_sgRNA-Read-Counts.xlsx`.

Unlike the Liu re-analysis, these are **raw integer read counts** — the per-library
sgRNA counts the authors produced with `mageck count` and then normalised, not a
published normalised table. Both methods in the comparison start from them, so
the head-to-head is a comparison of models on identical input.

## What the script derives

The workbook holds four sheets, one per library set. Every sheet name is prefixed
`CRISPRa`, but Calabrese is the CRISPRa library and Dolcetto the CRISPRi library
(Methods, and the GEO sample titles), so modality is taken from the library:

| Sheet | Prepared as | Guides | Reads |
|---|---|---|---|
| `CRISPRa.CalabreseSetA.count` | `CRISPRa_SetA` | 56,762 | 54.7M |
| `CRISPRa.CalabreseSetB.count` | `CRISPRa_SetB` | 56,476 | 47.5M |
| `CRISPRa.DolcettoSetA.count` | `CRISPRi_SetA` | 57,050 | 61.1M |
| `CRISPRa.DolcettoSetB.count` | `CRISPRi_SetB` | 57,011 | 52.8M |

Each sheet carries a plasmid library plus twelve sorted libraries — two donors ×
{IL-2, IFN-γ} × {low, high, unsorted}. The script drops the plasmid column, which
is the pre-transduction pool and plays no part in the high-versus-low contrast,
decodes donor, cytokine and bin from the column names, and records each library's
column sum as its beta-binomial denominator. No normalisation is applied.

Set A and Set B target the same ~18,800 genes with three guides each and share no
sgRNA sequence at all — including the 496 nontargeting sgRNAs, which are
different sequences in the two sets (Set A of Calabrese and Set A of Dolcetto do
share theirs). A gene therefore has six guides, three from each separately
sequenced pool, which is what the published pipeline's merge step is for and what
BARCS instead handles at the gene level.
