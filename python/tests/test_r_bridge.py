"""Unit checks on the rpy2 marshalling in `barcs`.

`test_matches_rscript.py` runs the whole pipeline both ways and compares the
numbers. This file pins the individual conversions that can break quietly,
and covers the functions that exist only because `barcs` wraps R rather than
reimplementing it.

There is no statistical comparison to make -- the wrapper *is* R, so a
disagreement would be a conversion bug, not a numerical one. What can
plausibly go wrong:

*   a guide-by-sample matrix arriving transposed, because R fills matrices
    column-major;
*   a pandas Categorical losing its level order, which would silently change
    the treatment coding and hence what `term=` refers to;
*   R's `NA` and logical columns coming back as something other than NaN and
    bool;
*   attributes such as the control scale, which do not survive the data.frame
    conversion and have to be fetched separately.

Skipped when rpy2 or R is unavailable.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

pytest.importorskip("rpy2", reason="rpy2 is not installed")

from importlib import import_module  # noqa: E402

import barcs as barcs_r  # noqa: E402

# `barcs.bbreg` the attribute is the re-exported *function*, which shadows the
# submodule of the same name, so the module has to be imported explicitly to
# reach its private conversion helpers.
barcs_module = import_module("barcs.bbreg")
from test_matches_rscript import build_fixture  # noqa: E402


@pytest.fixture(scope="module")
def fixture(tmp_path_factory):
    return build_fixture(tmp_path_factory.mktemp("bridge_fixture"))


@pytest.fixture(scope="module")
def screen(fixture) -> pd.DataFrame:
    return barcs_r.bb_screen(
        fixture["counts"],
        fixture["samples"],
        "~ dose + batch",
        term="dose",
        totals=fixture["totals"],
        guide=fixture["guide"],
        gene=fixture["gene"],
    )


def test_session_sources_the_r_implementation():
    import rpy2.robjects as ro

    available = set(ro.r("ls()"))
    # Everything the module claims to expose must actually be there, including
    # the functions the pure-Python port does not implement.
    for name in (
        "bbreg",
        "bb_contrast",
        "bb_screen",
        "bb_calibrate_controls",
        "bb_moderate_dispersion",
        "bb_gene_original",
        "bb_gene_normal",
        "bb_gene_partial_pool",
        "bb_gene_eb_moderate",
    ):
        assert name in available, f"{name} was not sourced from R/bbreg.R"


def test_single_fit_returns_a_usable_table(fixture):
    fit = barcs_r.bbreg(
        fixture["counts"][0], fixture["totals"], "~ dose + batch", fixture["samples"]
    )
    table = fit.coefficient_table

    # Coefficient names are R's, which is what `term=` is matched against.
    assert list(table.index) == ["(Intercept)", "dose", "batchb2", "batchb3"]
    assert list(table.columns) == ["estimate", "std_error", "t_value", "df", "p_value"]
    assert np.isfinite(table["estimate"]).all()

    contrast = barcs_r.bb_contrast(fit, {"dose": 2.0})
    # Two dose steps is exactly twice one, which also confirms the named
    # contrast reached R with its names attached rather than as a bare vector.
    assert contrast["estimate"][0] == pytest.approx(2 * table.loc["dose", "estimate"])


def test_counts_matrix_is_not_transposed(fixture, screen):
    """A transposed matrix is the failure this bridge is most likely to have.

    It would not raise: with 96 guides and 15 samples R would error on the
    shape, but on a square-ish input it would quietly fit the wrong thing. The
    row count is checked, and then that the values correspond to the guides
    they claim to.
    """
    assert len(screen) == fixture["counts"].shape[0]
    assert list(screen["guide"]) == list(fixture["guide"])

    # mean_cpm is computed inside R from the matrix it received, so comparing
    # it against the same quantity computed here in Python from the array we
    # sent proves the orientation survived the trip.
    expected = (fixture["counts"] / fixture["totals"] * 1e6).mean(axis=1)
    np.testing.assert_allclose(screen["mean_cpm"].to_numpy(), expected, rtol=1e-12)


def test_factor_levels_survive_conversion(fixture):
    """Treatment coding must match a native R session, or `term=` shifts."""
    import rpy2.robjects as ro

    converted = barcs_module._to_r(fixture["samples"])
    levels = list(ro.r["levels"](converted.rx2("batch")))
    assert levels == ["b1", "b2", "b3"], "factor level order changed crossing into R"


def test_types_come_back_as_python_types(screen):
    assert screen["converged"].dtype == bool
    assert screen["estimate"].dtype == float
    # BH-adjusted FDR is a probability; a conversion that mangled it would
    # most likely produce something outside the unit interval.
    finite = screen["fdr"].dropna()
    assert ((finite >= 0) & (finite <= 1)).all()


def test_control_calibration_returns_its_attribute(fixture, screen):
    control = np.asarray(fixture["gene"]) == "control"
    calibrated = barcs_r.bb_calibrate_controls(screen, control=control)

    scale = calibrated.attrs["control_scale"]
    assert scale is not None, "control_scale attribute was lost in conversion"
    assert scale >= 1.0, "min_scale defaults to 1, so calibration cannot be liberalising"
    # Calibration divides the t statistics by the scale.
    np.testing.assert_allclose(
        calibrated["t_value"].to_numpy(),
        screen["t_value"].to_numpy() / scale,
        rtol=1e-12,
    )


@pytest.mark.parametrize(
    "name",
    [
        "bb_gene_original",
        "bb_gene_normal",
        "bb_gene_partial_pool",
        "bb_gene_eb_moderate",
    ],
)
def test_every_gene_statistic_is_reachable(screen, name):
    """The four gene statistics, including the three the port does not have.

    This is the reason the bridge exists, so each one is actually called rather
    than assumed to work.
    """
    function = getattr(barcs_r, name)
    genes = function(screen)

    assert len(genes) > 0
    for column in ("gene", "estimate", "p_value", "fdr"):
        assert column in genes.columns, f"{name} did not return a `{column}` column"
    probabilities = genes["p_value"].dropna()
    assert ((probabilities >= 0) & (probabilities <= 1)).all()


def test_dispersion_moderation_is_reachable(screen):
    """Also absent from the pure-Python module."""
    moderated = barcs_r.bb_moderate_dispersion(screen)

    assert len(moderated) == len(screen)
    # Moderation keeps the unmoderated columns under a prefix, which is the
    # clearest signal the R function really ran rather than the frame being
    # passed through.
    assert any(column.startswith("unmoderated_") for column in moderated.columns)
    assert moderated.attrs["prior_df"] is not None
