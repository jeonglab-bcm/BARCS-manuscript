"""BARCS: beta-binomial analysis and regression for CRISPR screens.

A Python port of the guide-level regression layer in ``R/bbreg.R``. See that
module's docstring for what is and is not covered, and
``python/README.md`` for how the two implementations are kept in step.

    from barcs import bbreg, bb_screen

    fit = bbreg(counts_for_one_guide, library_totals, "~ dose + batch", samples)
    print(fit.summary())
"""

from barcs.bbreg import (
    BBRegResult,
    bb_calibrate_controls,
    bb_contrast,
    bb_gene_original,
    bb_screen,
    bbreg,
    benjamini_hochberg,
)

__all__ = [
    "BBRegResult",
    "bbreg",
    "bb_contrast",
    "bb_screen",
    "bb_calibrate_controls",
    "bb_gene_original",
    "benjamini_hochberg",
]

__version__ = "0.1.0"
