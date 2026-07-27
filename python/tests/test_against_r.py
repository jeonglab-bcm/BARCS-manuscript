"""Check the Python port against the R implementation it was ported from.

The port only earns trust by agreeing with the original, so this builds a
fixture screen, fits it both ways, and requires the two to match.

On the tolerance
----------------
Coefficient estimates agree to about 1e-12, which is double-precision noise.
Standard errors agree to about 4e-10, and the reason is worth stating because
it is not a defect in either implementation.

Both estimate ``rho`` by solving ``pearson(rho) = df_residual``. R uses
``uniroot(tol = 1e-10)``, which on the fixture stops at a point where the
Pearson statistic is 11.0000000083 against a target of 11; SciPy's ``brentq``,
asked for a tighter tolerance, reaches 10.999999999999776. So the two land on
``rho`` values that differ by 7.5e-10 relative -- and the Python one is the
more accurate root. That difference flows into the working weights and hence
into the standard errors.

That 7.5e-10 in ``rho`` is the seed of every other difference, and it grows a
little at each step that divides one noisy quantity by another. Standard
errors land at 4e-10, the negative-control scale -- a ratio of quantiles of
``estimate / std_error`` -- at 1e-8, and the calibrated FDR at 3e-8.

The threshold is therefore 1e-7: still eight significant figures of agreement
on every reported quantity, and honest about the fact that two root finders
with different stopping rules are not expected to agree to the last bit.
Tightening it further would mean deliberately making the Python solver worse
in order to imitate R's stopping rule.

The fixture is deliberately awkward, because agreement on easy data proves
little. It includes a continuous covariate crossed with a three-level factor,
guides spanning three orders of magnitude of abundance, genuinely
overdispersed counts, and a set of true nulls to calibrate against.

Skipped when Rscript or the CB2 environment is unavailable, so the pure-Python
suite still runs anywhere:

    pixi run test-python           # Python only
    pixi run test-python-vs-r      # includes this file
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

from barcs import (  # noqa: E402
    bb_calibrate_controls,
    bb_contrast,
    bb_gene_original,
    bb_screen,
    bbreg,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
TOLERANCE = 1e-7  # see the module docstring: R's uniroot stopping rule sets this

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
    # and scarce guides. The floor is chosen empirically for this seed as the
    # widest range in which every guide still converges; the assertion in
    # test_screen_matches_r fails loudly if a future edit breaks that.
    #
    # It matters because an IRLS that stopped at `maxit` is not obliged to stop
    # in the same place in two languages -- on an earlier fixture the two
    # returned +66 and +49 for the same unconverged guide -- and one such guide
    # contaminates every screen-wide quantity: the BH correction, the control
    # scale, and the gene aggregation are all computed over all guides at once.
    #
    # Convergence is not simply a function of depth. A guide with 448 reads can
    # fail where one with 237 succeeds, because what fails is the fixed-point
    # iteration rather than the arithmetic, so this floor is a property of this
    # seed and not a general threshold. Sparse-guide behaviour is covered by
    # test_sparse_guides_are_reported_not_raised instead.
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
        "design matrix columns differ between formulaic and R's model.matrix"
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

    The two implementations part company here, and it is recorded rather than
    smoothed over. On a guide with only a handful of reads spread across many
    samples, R's `bbreg` errors and `bb_screen` catches it into an all-NA row;
    the Python IRLS sometimes still returns a value. Neither answer is
    meaningful at that abundance.

    What both guarantee, and what this pins down, is narrower than it might
    seem: one unusable guide never aborts the screen, and a guide with no reads
    at all is always reported as unusable rather than as a result.

    A guide with a *handful* of reads is not covered by that. Python often
    converges on one and returns a finite estimate where R declines to fit it,
    so `converged` being True on such a guide means only that the iteration
    stopped moving -- not that the answer means anything. Filter on abundance,
    via `min_total_count` or afterwards; do not rely on `converged` to do it,
    and do not compare the two implementations on guides like these.
    """
    rng = np.random.default_rng(11)
    samples = pd.DataFrame({"dose": np.tile([-1.0, 0.0, 1.0], 4)})
    totals = np.full(12, 50_000.0)

    counts = rng.binomial(50_000, 0.0004, size=(4, 12)).astype(float)
    counts[2] = 0.0  # never observed
    counts[3] = np.array([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], dtype=float)  # two reads

    screen = bb_screen(counts, samples, "~ dose", term="dose", totals=totals)

    assert len(screen) == 4, "a sparse guide must not drop rows from the screen"

    # The all-zero guide falls below `min_total_count` and must never be
    # presented as an estimate.
    assert not bool(screen["converged"][2])
    assert not np.isfinite(screen["estimate"][2])

    # The two-read guide is the documented grey area: whatever it reports, it
    # must not have derailed the screen, and its mean CPM must still be
    # recorded so it can be filtered on abundance afterwards.
    assert np.isfinite(screen["mean_cpm"][3])

    # The ordinary guides are unaffected by their sparse neighbours.
    assert bool(screen["converged"][0])
    assert np.isfinite(screen["estimate"][0])
