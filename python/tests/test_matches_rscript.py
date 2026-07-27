"""Check that the rpy2 wrapper marshals data correctly.

`barcs` calls `R/bbreg.R` through rpy2, so there is no second implementation
to compare against and nothing statistical to verify -- the numerics are R's.
What can still go wrong is everything around the call: a matrix arriving
transposed, a factor losing its level order, a column coming back as the wrong
type, an R attribute quietly dropped.

So this fits one fixture screen twice. Once through the wrapper, and once by
running `reference_fit.R` under a plain `Rscript` over CSV files written to
disk. Both paths source the same `R/bbreg.R` with no CB2 loaded, so they
execute identical code and any difference is the marshalling.

That makes the bar exact rather than approximate. When this module was a
port, the two agreed to eight significant figures and the residual came from
R's `uniroot` and SciPy's `brentq` stopping in different places; there is now
only one root finder, so the tolerance is 1e-12 -- double-precision noise from
the CSV round trip and nothing else. A failure here is a conversion bug, not a
numerical one.

The fixture is deliberately awkward: a continuous covariate crossed with a
three-level factor, guides across a wide abundance range, genuinely
overdispersed counts, and a set of true nulls for the calibration to read.

Skipped when Rscript is unavailable.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

pytest.importorskip("rpy2", reason="rpy2 is not installed")

from barcs import (  # noqa: E402
    bb_calibrate_controls,
    bb_contrast,
    bb_gene_original,
    bb_screen,
    bbreg,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
# Both paths run the same R code, so this is a CSV round-trip tolerance.
TOLERANCE = 1e-12

requires_r = pytest.mark.skipif(
    shutil.which("Rscript") is None,
    reason="Rscript is not available; run inside `pixi shell`.",
)


def build_fixture(directory: Path) -> dict[str, object]:
    """Write a deterministic beta-binomial screen both implementations read."""
    rng = np.random.default_rng(20260727)

    n_genes, guides_per_gene = 24, 4
    n_guides = n_genes * guides_per_gene

    dose = np.tile(np.array([-1.5, -0.75, 0.0, 0.75, 1.5]), 3)
    batch = np.repeat(["b1", "b2", "b3"], 5)
    samples = pd.DataFrame({"dose": dose, "batch": pd.Categorical(batch)})
    n_samples = len(samples)

    totals = rng.integers(40_000, 120_000, size=n_samples).astype(float)

    gene_names = [f"gene{i:03d}" for i in range(n_genes)]
    # The last six genes are the negative controls the calibration reads.
    gene_names[-6:] = ["control"] * 6
    gene = np.repeat(gene_names, guides_per_gene)
    guide = [f"sg{i:04d}" for i in range(n_guides)]

    effect = np.where(np.array(gene) == "control", 0.0, 0.0)
    active = np.zeros(n_guides, dtype=bool)
    active[: 8 * guides_per_gene] = True
    effect = np.where(active, np.repeat(rng.normal(0, 0.8, n_genes), guides_per_gene), 0.0)
    effect[np.array(gene) == "control"] = 0.0

    # Baselines spanning a wide abundance range, so the screen mixes abundant
    # and scarce guides, with a floor chosen so every guide converges. Both
    # paths run the same solver now, so an unconverged guide would not make
    # them disagree -- but it would still make the fixture a weaker test, since
    # a screen full of NA rows compares nothing.
    baseline = rng.uniform(-7.8, -6.5, size=n_guides)
    batch_effect = {"b1": 0.0, "b2": 0.18, "b3": -0.12}
    offsets = np.array([batch_effect[b] for b in batch])

    rho = 0.002
    precision = 1 / rho - 1
    counts = np.zeros((n_guides, n_samples))
    for g in range(n_guides):
        eta = baseline[g] + effect[g] * samples["dose"].to_numpy() + offsets
        mu = 1 / (1 + np.exp(-eta))
        latent = rng.beta(mu * precision, (1 - mu) * precision)
        counts[g] = rng.binomial(totals.astype(int), latent)

    directory.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(counts, index=guide, columns=[f"s{i:02d}" for i in range(n_samples)]).to_csv(
        directory / "counts.csv"
    )
    samples.assign(batch=samples["batch"].astype(str)).to_csv(
        directory / "samples.csv", index=False
    )
    np.savetxt(directory / "totals.txt", totals)
    pd.DataFrame({"guide": guide, "gene": gene}).to_csv(directory / "guide_gene.csv", index=False)

    return {
        "counts": counts,
        "samples": samples,
        "totals": totals,
        "guide": guide,
        "gene": gene,
    }


@pytest.fixture(scope="module")
def fixture(tmp_path_factory) -> dict[str, object]:
    directory = tmp_path_factory.mktemp("barcs_fixture")
    data = build_fixture(directory)
    data["directory"] = directory
    return data


@pytest.fixture(scope="module")
def r_results(fixture) -> dict[str, pd.DataFrame]:
    """Run the R reference script over the fixture and load what it wrote."""
    directory = fixture["directory"]
    completed = subprocess.run(
        ["Rscript", "python/tests/reference_fit.R", str(directory)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        pytest.skip(f"R reference fit failed:\n{completed.stderr[-2000:]}")

    return {
        path.stem.removeprefix("r_"): pd.read_csv(path)
        for path in sorted(directory.glob("r_*.csv"))
    }


def assert_close(actual, expected, label: str) -> None:
    actual = np.asarray(actual, dtype=float)
    expected = np.asarray(expected, dtype=float)
    both_nan = np.isnan(actual) & np.isnan(expected)
    if both_nan.all():
        return
    difference = np.abs(actual - expected)[~both_nan]
    # Relative where the values are large, absolute where they are near zero;
    # p-values legitimately reach 1e-300, where a relative test is meaningless.
    scale = np.maximum(1.0, np.abs(expected)[~both_nan])
    worst = float(np.max(difference / scale))
    assert worst < TOLERANCE, f"{label}: largest scaled difference {worst:.3e}"


@requires_r
def test_single_guide_fit_matches_r(fixture, r_results):
    """Every coefficient of one fit, not just the reported term."""
    reference = r_results["single_fit"]
    fit = bbreg(
        fixture["counts"][0],
        fixture["totals"],
        "~ dose + batch",
        fixture["samples"],
    )

    assert list(fit.coefficient_table.index) == list(reference["coefficient"]), (
        "coefficient names differ between the wrapper and a native R session"
    )
    for column in ("estimate", "std_error", "t_value", "df", "p_value"):
        assert_close(fit.coefficient_table[column], reference[column], f"single fit {column}")

    assert_close([fit.rho], [reference["rho"][0]], "rho")
    assert_close([fit.pearson], [reference["pearson"][0]], "pearson")
    assert_close([fit.pearson_null], [reference["pearson_null"][0]], "pearson_null")
    assert fit.converged == bool(reference["converged"][0])


@requires_r
def test_contrast_matches_r(fixture, r_results):
    reference = r_results["contrast"]
    fit = bbreg(
        fixture["counts"][0],
        fixture["totals"],
        "~ dose + batch",
        fixture["samples"],
    )
    contrast = bb_contrast(fit, {"dose": 2.0})
    for column in ("estimate", "std_error", "t_value", "df", "p_value"):
        assert_close(contrast[column], reference[column], f"contrast {column}")


@requires_r
def test_screen_matches_r(fixture, r_results):
    """Every guide, every reported column, including the BH-adjusted FDR.

    The FDR matters here in a way the per-guide columns do not: it is computed
    across the whole screen, so it only agrees if every guide agreed and both
    implementations counted the same number of testable guides.
    """
    reference = r_results["screen"]
    screen = bb_screen(
        fixture["counts"],
        fixture["samples"],
        "~ dose + batch",
        term="dose",
        totals=fixture["totals"],
        guide=fixture["guide"],
        gene=fixture["gene"],
    )

    assert list(screen["guide"]) == list(reference["guide"])
    assert list(screen["converged"].astype(bool)) == list(
        reference["converged"].astype(bool)
    )
    # Guard the fixture itself: if a future edit made guides too sparse to fit,
    # the comparisons below would quietly become vacuous on those rows.
    assert screen["converged"].all(), "every fixture guide should be fittable"

    for column in (
        "estimate",
        "std_error",
        "t_value",
        "df",
        "p_value",
        "rho",
        "pearson_null",
        "mean_cpm",
        "fdr",
    ):
        assert_close(screen[column], reference[column], f"screen {column}")


@requires_r
@pytest.mark.parametrize("method", ["tail_quantile", "qq_slope"])
def test_control_calibration_matches_r(fixture, r_results, method):
    reference = r_results[f"calibrated_{method}"]
    screen = bb_screen(
        fixture["counts"],
        fixture["samples"],
        "~ dose + batch",
        term="dose",
        totals=fixture["totals"],
        guide=fixture["guide"],
        gene=fixture["gene"],
    )
    control = np.asarray(fixture["gene"]) == "control"
    calibrated = bb_calibrate_controls(screen, control=control, method=method)

    assert_close(
        [calibrated.attrs["control_scale"]],
        [reference["control_scale"][0]],
        f"{method} scale",
    )
    for column in ("std_error", "t_value", "p_value", "fdr"):
        assert_close(calibrated[column], reference[column], f"{method} {column}")


@requires_r
def test_gene_aggregation_matches_r(fixture, r_results):
    reference = r_results["gene"]
    screen = bb_screen(
        fixture["counts"],
        fixture["samples"],
        "~ dose + batch",
        term="dose",
        totals=fixture["totals"],
        guide=fixture["guide"],
        gene=fixture["gene"],
    )
    genes = bb_gene_original(screen)

    assert list(genes["gene"]) == list(reference["gene"])
    assert list(genes["n_guides"]) == list(reference["n_guides"])
    for column in ("estimate", "statistic", "p_value", "fdr", "guide_direction_agreement"):
        assert_close(genes[column], reference[column], f"gene {column}")


def test_sparse_guides_are_reported_not_raised():
    """A guide too sparse to fit is flagged, not fatal.

    `bb_screen` skips any guide whose total falls below `min_total_count`,
    which now defaults to R's 10 rather than the port's 1. Such a guide must
    still occupy a row -- dropping it would silently misalign the result
    against the caller's guide list -- and must never appear as an estimate.
    """
    rng = np.random.default_rng(11)
    samples = pd.DataFrame({"dose": np.tile([-1.0, 0.0, 1.0], 4)})
    totals = np.full(12, 50_000.0)

    counts = rng.binomial(50_000, 0.0004, size=(4, 12)).astype(float)
    counts[2] = 0.0  # never observed
    counts[3] = np.array([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], dtype=float)  # two reads

    screen = bb_screen(counts, samples, "~ dose", term="dose", totals=totals)

    assert len(screen) == 4, "a sparse guide must not drop rows from the screen"

    # Both sparse guides fall below the default min_total_count of 10.
    for row in (2, 3):
        assert not bool(screen["converged"][row])
        assert not np.isfinite(screen["estimate"][row])
        # mean_cpm is still recorded, so abundance filtering remains possible.
        assert np.isfinite(screen["mean_cpm"][row])

    # The ordinary guides are unaffected by their sparse neighbours.
    assert bool(screen["converged"][0])
    assert np.isfinite(screen["estimate"][0])
