# Chronos/Tzelepis benchmark sources

The two compressed text files were extracted from Data S2 (`mmc6.zip`) of
Tzelepis et al.:

- Europe PMC supplementary bundle:
  `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5081405/supplementaryFiles`
- `Tzelepis_HT29_counts.txt.gz` is the gzip-compressed
  `grna_counts_HT29C-clone3_d7-25-HT1080.txt`.
- `Tzelepis_HT29_expression_FPKM.txt.gz` is the gzip-compressed
  `RNA-seq_5AMLs-HT1080-HT29_FPKM_Kallisto.txt`.

The HDF5 effect matrices and control lists are individual files from the
Chronos Figshare article `14067047`:

| Local file | Figshare file ID |
|---|---:|
| `GeneFitnessEffect_ChronosJoint_Tzelepis.hdf5` | 26548532 |
| `GeneFitnessEffect_MAGeCK_Tzelepis.hdf5` | 26548541 |
| `GeneFitnessEffect_BAGEL2_Tzelepis.hdf5` | 30847513 |
| `ReferenceEssentials.csv` | 26548550 |
| `ReferenceNonEssentials.csv` | 26548553 |
| `DepMapSampleInfo20Q2.csv` | 26548499 |

An individual Figshare file can be downloaded from
`https://ndownloader.figshare.com/files/FILE_ID`.

`data/derived/HT29_DepMap20Q2_CNV.tsv` contains the ACH-000552 column extracted
from Figshare file 26548484 (`CCLEGeneCopyNumber20Q2.hdf5`). The source matrix
stores gene names in `dim_1`, DepMap IDs in `dim_0`, and values in `data`.
Likewise, `data/derived/HT29_DepMap20Q2_expression.tsv` contains ACH-000552
expression extracted from Figshare file 26548487
(`CCLEGeneExpression20Q2.hdf5`). This is the expression source used by the
Chronos analysis notebook to define genes below 0.5 expression.
