# Manuscript-wide method palette.
#
# A method keeps the same base colour in every figure. Variants such as raw
# versus calibrated or uncorrected versus CNV-corrected reuse the base colour
# and are distinguished by line type, point shape, or fill.
barcs_method_colours <- c(
  BARCS = "#0072B2",
  `BARCS-original` = "#0072B2",
  `BARCS-ST` = "#0072B2",
  `BARCS-MOD` = "#0072B2",
  `BARCS-moderated` = "#0072B2",
  `BARCS-NORM` = "#7A3E9D",
  `BARCS-partial` = "#009E73",
  `BARCS-RE` = "#009E73",
  `BARCS-EB` = "#D55E00",
  `BARCS-RE-EB` = "#D55E00",
  MAGeCK = "#D55E00",
  `MAGeCK-RRA` = "#D55E00",
  `MAGeCK-MLE` = "#D55E00",
  `Liang RRA` = "#CC79A7",
  edgeR = "#009E73",
  DESeq2 = "#CC79A7",
  `limma-voom` = "#E69F00",
  Chronos = "#009E73",
  `Published MAGeCK` = "#CC79A7",
  BAGEL2 = "#E69F00",
  Waterbear = "#7A3E9D",
  MAUDE = "#A33D3D",
  `Naive binomial` = "#767676",
  CB2 = "#56B4E9",
  Validation = "#333333",
  # Okabe-Ito is exhausted by the entries above, so CRISPhieRmix reuses the
  # BARCS-NORM/Waterbear violet. It never appears in a figure with either of
  # them; if that changes, one of the three needs a new hue.
  CRISPhieRmix = "#7A3E9D"
)
