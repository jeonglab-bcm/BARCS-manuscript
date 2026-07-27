"""BARCS: beta-binomial analysis and regression for CRISPR screens.

A thin rpy2 wrapper around the guide-level regression layer in
``R/bbreg.R``. See that module's docstring for why this calls into R rather
than reimplementing it, and ``python/README.md`` for setup and usage.

    from barcs import bbreg, bb_screen

    fit = bbreg(counts_for_one_guide, library_totals, "~ dose + batch", samples)
    print(fit.summary())
"""

from barcs.bbreg import (
    BBRegResult,
    bb_calibrate_controls,
    bb_contrast,
    bb_gene_consistency,
    bb_gene_eb_moderate,
    bb_gene_normal,
    bb_gene_original,
    bb_gene_partial_pool,
    bb_moderate_dispersion,
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
    "bb_moderate_dispersion",
    "bb_gene_original",
    "bb_gene_normal",
    "bb_gene_consistency",
    "bb_gene_partial_pool",
    "bb_gene_eb_moderate",
    "benjamini_hochberg",
]

__version__ = "0.2.0"

