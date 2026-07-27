"""Beta-binomial regression with t-based coefficient tests.

A Python port of the regression layer in ``R/bbreg.R``. It fits the mean model

.. math::

    \\operatorname{logit}(\\mu_i) = x_i' \\beta

under the beta-binomial variance

.. math::

    \\operatorname{Var}(K_i) = n_i \\mu_i (1 - \\mu_i)
                               \\{1 + (n_i - 1) \\rho\\}.

Coefficients come from feasible IRLS, ``rho`` from the Pearson estimating
equation, and coefficient tests use a Student t reference on the residual
degrees of freedom.

Why a port rather than a wrapper
--------------------------------
Calling R from Python would have kept one implementation, at the cost of
making every Python user install R and the package. This is instead a
standalone implementation with no R dependency at all -- but that only means
something if it agrees with the original, so ``python/tests/`` checks the two
against each other on shared fixtures and requires agreement to 1e-10.
``pixi run test-python-vs-r`` runs that comparison.

Read the R file for the statistical reasoning; it is the reference text and is
not repeated here. Comments below cover what is specific to reproducing it in
Python, which is mostly places where a NumPy or SciPy default differs from
R's.

What is ported
--------------
The guide-level inference path, end to end: :func:`bbreg`,
:func:`bb_contrast`, :func:`bb_screen`, :func:`bb_calibrate_controls`, and the
``original`` guide-to-gene statistic :func:`bb_gene_original`.

Not ported: the exchangeable-normal, partial-pooling, and empirical-Bayes
guide-to-gene statistics, and :func:`bb_moderate_dispersion`. Those are what
the manuscript's later comparisons turn on, and every committed result
involving them came from the R implementation. Use R for those rather than
assuming this module covers them.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Sequence

import numpy as np
import pandas as pd
from scipy import optimize, stats
from scipy.special import expit

__all__ = [
    "BBRegResult",
    "bbreg",
    "bb_contrast",
    "bb_screen",
    "bb_calibrate_controls",
    "bb_gene_original",
    "benjamini_hochberg",
]

# R's `.Machine$double.eps` and `double.xmin`.
_DOUBLE_EPS = np.finfo(float).eps
_DOUBLE_XMIN = np.finfo(float).tiny


# ---------------------------------------------------------------------------
# Small numerical helpers where R and NumPy defaults differ
# ---------------------------------------------------------------------------


def benjamini_hochberg(p_value: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values, matching ``p.adjust(method="BH")``.

    Written out rather than taken from statsmodels for two reasons: it avoids a
    dependency for twelve lines of arithmetic, and R's version propagates NA
    while ranking against ``n`` = the number of *non-missing* p-values. A guide
    that failed to converge must not shrink the multiple-testing correction for
    the ones that did.
    """
    p_value = np.asarray(p_value, dtype=float)
    out = np.full(p_value.shape, np.nan)
    finite = np.isfinite(p_value)
    n = int(finite.sum())
    if n == 0:
        return out

    values = p_value[finite]
    order = np.argsort(-values, kind="stable")  # decreasing, as R does
    ranks = np.arange(n, 0, -1, dtype=float)
    # Cumulative minimum of (n / i) * p over the decreasing order, clipped at 1.
    adjusted = np.minimum.accumulate(np.minimum(n / ranks * values[order], 1.0))
    restored = np.empty(n)
    restored[order] = adjusted
    out[finite] = restored
    return out


def _quantile_type8(x: np.ndarray, probs: np.ndarray) -> np.ndarray:
    """R's ``quantile(type = 8)``.

    NumPy's default is R's type 7. The control calibration in the R code asks
    for type 8 explicitly -- the median-unbiased definition -- and the two give
    visibly different answers in the far tail, which is exactly where the
    empirical null scale is read off. NumPy spells it ``median_unbiased``.
    """
    return np.quantile(np.asarray(x, dtype=float), probs, method="median_unbiased")


# ---------------------------------------------------------------------------
# Validation and design construction
# ---------------------------------------------------------------------------


def _validate_response(count: np.ndarray, total: np.ndarray) -> None:
    if count.shape != total.shape or count.size < 2:
        raise ValueError("`count` and `total` must have the same length (at least two).")
    if not (np.isfinite(count).all() and np.isfinite(total).all()):
        raise ValueError("`count` and `total` cannot contain missing or non-finite values.")
    if (total <= 0).any() or (count < 0).any() or (count > total).any():
        raise ValueError("Each observation must satisfy 0 <= `count` <= `total`, with `total` > 0.")
    tol = np.sqrt(_DOUBLE_EPS)
    if (np.abs(count - np.round(count)) >= tol).any() or (
        np.abs(total - np.round(total)) >= tol
    ).any():
        raise ValueError("`count` and `total` must contain integer-valued counts.")


def _r_style_name(column: str) -> str:
    """Rename a formulaic column to what R's ``model.matrix`` would call it.

    The two libraries build the same columns and disagree only on the labels:
    formulaic writes ``Intercept`` and ``batch[T.b2]`` where R writes
    ``(Intercept)`` and ``batchb2``. The labels are not cosmetic here --
    :func:`bb_screen` selects a coefficient by name through its ``term``
    argument, so leaving them alone would mean the same analysis needed a
    different ``term`` in each language, and every comparison against an R
    result would have to translate.
    """
    if column == "Intercept":
        return "(Intercept)"
    # `batch[T.b2]` -> `batchb2`; leaves plain numeric terms such as `dose`
    # untouched, and handles interactions by rewriting each `:`-joined part.
    return ":".join(
        re.sub(r"^(.*)\[T\.(.*)\]$", r"\1\2", part) for part in column.split(":")
    )


def _make_design(formula: str | np.ndarray, data: pd.DataFrame, n: int):
    """Build the model matrix, matching R's ``model.matrix`` conventions.

    A design matrix can be passed directly, which keeps this usable without
    ``formulaic``. When a formula string is given, formulaic parses it with
    treatment coding and the first factor level as the reference -- R's
    defaults -- so ``~ dose + batch`` produces the same columns in the same
    order as R.
    """
    if isinstance(formula, str):
        try:
            from formulaic import model_matrix
        except ImportError as exc:  # pragma: no cover - environment guard
            raise ImportError(
                "Formula input needs `formulaic`. Install it, or pass an "
                "already-built design matrix instead."
            ) from exc
        if not isinstance(data, pd.DataFrame) or len(data) != n:
            raise ValueError("`data` must be a data frame with one row per observation.")
        matrix = model_matrix(formula, data)
        names = [_r_style_name(column) for column in matrix.columns]
        x = np.asarray(matrix, dtype=float)
    else:
        x = np.asarray(formula, dtype=float)
        if x.ndim != 2 or x.shape[0] != n:
            raise ValueError("A design matrix must have one row per observation.")
        names = [f"x{i}" for i in range(x.shape[1])]

    if np.linalg.matrix_rank(x) < x.shape[1]:
        raise ValueError("The design matrix is not full rank; remove aliased covariates.")
    if x.shape[0] <= x.shape[1]:
        raise ValueError("The model needs more observations than fitted coefficients.")
    return x, names


# ---------------------------------------------------------------------------
# Dispersion
# ---------------------------------------------------------------------------


@dataclass
class _RhoFit:
    rho: float
    scale: float
    pearson: float
    pearson_null: float
    boundary: bool


def _estimate_rho(
    count: np.ndarray,
    total: np.ndarray,
    mu: np.ndarray,
    df_residual: int,
    upper: float = 1 - 1e-8,
) -> _RhoFit:
    binomial_variance = total * mu * (1 - mu)

    def pearson(rho: float) -> float:
        return float(
            np.sum((count - total * mu) ** 2 / (binomial_variance * (1 + (total - 1) * rho)))
        )

    q0 = pearson(0.0)
    # `pearson_null` is the Pearson chi-square under the binomial model.
    # Unlike `rho` it is not truncated at the lower boundary, which is what
    # makes it the right quantity for dispersion moderation to shrink.
    if not np.isfinite(q0):
        return _RhoFit(0.0, 1.0, q0, q0, True)
    if q0 <= df_residual:
        return _RhoFit(0.0, 1.0, q0, q0, False)

    qu = pearson(upper)
    if qu >= df_residual:
        return _RhoFit(upper, max(1.0, qu / df_residual), qu, q0, True)

    # R uses `uniroot(tol = 1e-10)`; brentq with the same bracket and a tighter
    # tolerance lands on the same root to well inside the comparison threshold.
    rho = float(
        optimize.brentq(lambda r: pearson(r) - df_residual, 0.0, upper, xtol=1e-12, rtol=8.9e-16)
    )
    return _RhoFit(rho, 1.0, pearson(rho), q0, False)


# ---------------------------------------------------------------------------
# Single-guide fit
# ---------------------------------------------------------------------------


@dataclass
class BBRegResult:
    """A fitted single-guide beta-binomial regression."""

    coefficients: pd.Series
    coefficient_table: pd.DataFrame
    covariance: np.ndarray
    fitted_values: np.ndarray
    linear_predictors: np.ndarray
    residuals: np.ndarray
    pearson: float
    pearson_null: float
    rho: float
    scale: float
    dispersion_boundary: bool
    df_residual: int
    rank: int
    count: np.ndarray = field(repr=False)
    total: np.ndarray = field(repr=False)
    design: np.ndarray = field(repr=False)
    weights: np.ndarray = field(repr=False)
    converged: bool = False
    iterations: int = 0

    def summary(self) -> str:
        lines = [
            "Beta-binomial regression",
            "",
            self.coefficient_table.to_string(),
            "",
            f"Beta-binomial intraclass correlation (rho): {self.rho:.6g}",
            f"Pearson statistic / residual df: {self.pearson:.6g} / {self.df_residual}",
            f"Converged: {self.converged} after {self.iterations} iterations",
        ]
        return "\n".join(lines)


def _initial_beta(x: np.ndarray, count: np.ndarray, total: np.ndarray, tolerance: float):
    """Binomial IRLS start, standing in for R's ``glm.fit``.

    R seeds the beta-binomial loop with an ordinary binomial GLM. Rather than
    depend on statsmodels for that one call, the equivalent IRLS is written out
    here. The starting value only has to be good, not identical -- the
    beta-binomial loop below iterates to a tolerance either way -- but keeping
    it the same algorithm avoids converging to a different local answer on
    badly behaved guides.
    """
    proportion = count / total
    mu = (count + 0.5) / (total + 1.0)
    beta = np.zeros(x.shape[1])
    for _ in range(50):
        eta = np.log(mu / (1 - mu))
        weight = total * mu * (1 - mu)
        working = eta + (proportion - mu) / (mu * (1 - mu))
        try:
            information = x.T @ (weight[:, None] * x)
            beta_new = np.linalg.solve(information, x.T @ (weight * working))
        except np.linalg.LinAlgError:
            return None
        if not np.isfinite(beta_new).all():
            return None
        change = np.max(np.abs(beta_new - beta) / np.maximum(1.0, np.abs(beta)))
        beta = beta_new
        eta = x @ beta
        mu = np.clip(expit(eta), 1e-10, 1 - 1e-10)
        if change < tolerance:
            break
    return beta


def bbreg(
    count: Sequence[float] | np.ndarray,
    total: Sequence[float] | np.ndarray,
    formula: str | np.ndarray,
    data: pd.DataFrame | None = None,
    maxit: int = 100,
    tolerance: float = 1e-8,
    mu_bound: float = 1e-8,
) -> BBRegResult:
    """Fit beta-binomial regression for one guide.

    Parameters mirror the R function. ``formula`` accepts either an R-style
    one-sided formula string such as ``"~ dose + batch"`` or a pre-built design
    matrix.
    """
    count = np.asarray(count, dtype=float)
    total = np.asarray(total, dtype=float)
    _validate_response(count, total)

    x, names = _make_design(formula, data, count.size)
    rank = x.shape[1]
    df_residual = x.shape[0] - rank

    beta = _initial_beta(x, count, total, tolerance)
    if beta is None or not np.isfinite(beta).all():
        pooled = (count.sum() + 0.5) / (total.sum() + 1.0)
        beta = np.zeros(rank)
        beta[0] = np.log(pooled / (1 - pooled))

    converged = False
    rho_fit = _RhoFit(0.0, 1.0, np.nan, np.nan, False)
    iteration = 0
    for iteration in range(1, maxit + 1):
        eta = x @ beta
        mu = np.clip(expit(eta), mu_bound, 1 - mu_bound)
        rho_fit = _estimate_rho(count, total, mu, df_residual)

        working_response = eta + (count / total - mu) / (mu * (1 - mu))
        working_weight = total * mu * (1 - mu) / (1 + (total - 1) * rho_fit.rho)

        try:
            information = x.T @ (working_weight[:, None] * x)
            beta_new = np.linalg.solve(information, x.T @ (working_weight * working_response))
        except np.linalg.LinAlgError as exc:
            raise ValueError(
                "The IRLS update was singular; inspect sparse counts and the design."
            ) from exc
        if not np.isfinite(beta_new).all():
            raise ValueError(
                "The IRLS update was singular; inspect sparse counts and the design."
            )

        change = np.max(np.abs(beta_new - beta) / np.maximum(1.0, np.abs(beta)))
        beta = beta_new
        if change < tolerance:
            converged = True
            break

    eta = x @ beta
    mu = np.clip(expit(eta), mu_bound, 1 - mu_bound)
    rho_fit = _estimate_rho(count, total, mu, df_residual)
    working_weight = total * mu * (1 - mu) / (1 + (total - 1) * rho_fit.rho)

    information = x.T @ (working_weight[:, None] * x)
    # R inverts via `chol2inv(chol(.))`. Cholesky here too: the information
    # matrix is positive definite by construction, and using the same
    # factorisation keeps the rounding identical rather than merely close.
    factor = np.linalg.cholesky(information)
    identity = np.eye(rank)
    unscaled = np.linalg.solve(factor.T, np.linalg.solve(factor, identity))
    covariance = rho_fit.scale * unscaled

    standard_error = np.sqrt(np.diag(covariance))
    statistic = beta / standard_error
    p_value = 2 * stats.t.cdf(-np.abs(statistic), df=df_residual)

    table = pd.DataFrame(
        {
            "estimate": beta,
            "std_error": standard_error,
            "t_value": statistic,
            "df": float(df_residual),
            "p_value": p_value,
        },
        index=names,
    )

    return BBRegResult(
        coefficients=pd.Series(beta, index=names),
        coefficient_table=table,
        covariance=covariance,
        fitted_values=mu,
        linear_predictors=eta,
        residuals=count / total - mu,
        pearson=rho_fit.pearson,
        pearson_null=rho_fit.pearson_null,
        rho=rho_fit.rho,
        scale=rho_fit.scale,
        dispersion_boundary=rho_fit.boundary,
        df_residual=df_residual,
        rank=rank,
        count=count,
        total=total,
        design=x,
        weights=working_weight,
        converged=converged,
        iterations=iteration,
    )


def bb_contrast(
    fit: BBRegResult,
    contrast: Sequence[float] | dict[str, float] | np.ndarray,
    null: float = 0.0,
) -> pd.DataFrame:
    """Test a linear contrast of the fitted coefficients."""
    names = list(fit.coefficients.index)
    if isinstance(contrast, dict):
        unknown = set(contrast) - set(names)
        if unknown:
            raise ValueError(f"A named contrast contains an unknown coefficient: {sorted(unknown)}")
        vector = np.zeros(len(names))
        for key, value in contrast.items():
            vector[names.index(key)] = value
    else:
        vector = np.asarray(contrast, dtype=float)
        if vector.size != len(names):
            raise ValueError("An unnamed contrast must have one value per coefficient.")
    if not np.isfinite(vector).all():
        raise ValueError("`contrast` must be numeric without missing values.")

    estimate = float(vector @ fit.coefficients.to_numpy())
    standard_error = float(np.sqrt(vector @ fit.covariance @ vector))
    statistic = (estimate - null) / standard_error
    return pd.DataFrame(
        {
            "estimate": [estimate],
            "std_error": [standard_error],
            "t_value": [statistic],
            "df": [float(fit.df_residual)],
            "p_value": [2 * stats.t.cdf(-abs(statistic), df=fit.df_residual)],
        }
    )


# ---------------------------------------------------------------------------
# Screen-level driver
# ---------------------------------------------------------------------------


def bb_screen(
    counts: np.ndarray | pd.DataFrame,
    data: pd.DataFrame,
    formula: str | np.ndarray,
    term: str,
    totals: Sequence[float] | None = None,
    guide: Sequence[str] | None = None,
    gene: Sequence[str] | None = None,
    min_total_count: int = 1,
    **kwargs: Any,
) -> pd.DataFrame:
    """Apply :func:`bbreg` guide by guide and report one coefficient.

    Returns one row per guide with the reported term's estimate, standard
    error, t statistic, degrees of freedom, p-value, and BH-adjusted FDR,
    alongside the guide's dispersion and mean CPM.
    """
    if isinstance(counts, pd.DataFrame):
        if guide is None:
            guide = list(counts.index)
        counts = counts.to_numpy(dtype=float)
    counts = np.asarray(counts, dtype=float)
    if counts.ndim != 2:
        raise ValueError("`counts` must be a guide-by-sample matrix.")
    n_guides, n_samples = counts.shape

    if len(data) != n_samples:
        raise ValueError("`data` must have one row per count-matrix column.")

    totals_array = counts.sum(axis=0) if totals is None else np.asarray(totals, dtype=float)
    if totals_array.size != n_samples:
        raise ValueError("`totals` must have one value per count-matrix column.")
    if (counts > totals_array[None, :]).any():
        raise ValueError("A guide count cannot exceed its sample's `total`.")
    # Checked here as well as per guide because the per-guide failure is
    # silent: every fit would return an all-NA row reading as a modelling
    # failure rather than a malformed argument. Size-factor normalization is
    # the usual way to arrive with non-integer totals.
    if (np.abs(totals_array - np.round(totals_array)) >= np.sqrt(_DOUBLE_EPS)).any():
        raise ValueError(
            "`totals` must be integer-valued library sizes; round them first. "
            "A beta-binomial denominator counts sequenced reads, so a fractional "
            "total has no likelihood."
        )

    if guide is None:
        guide = [f"guide_{i + 1}" for i in range(n_guides)]
    guide = list(guide)
    if len(guide) != n_guides or len(set(guide)) != n_guides:
        raise ValueError("`guide` must uniquely identify every row of `counts`.")
    if gene is not None and len(gene) != n_guides:
        raise ValueError("`gene` must have one value per guide.")

    _, names = _make_design(formula, data, n_samples)
    if term not in names:
        raise ValueError(f"`term` must be one model-matrix coefficient: {', '.join(names)}")

    records = []
    for i in range(n_guides):
        row = counts[i]
        mean_cpm = float(np.mean(row / totals_array * 1e6))
        blank = {
            "estimate": np.nan,
            "std_error": np.nan,
            "t_value": np.nan,
            "df": np.nan,
            "p_value": np.nan,
            "rho": np.nan,
            "pearson_null": np.nan,
            "mean_cpm": mean_cpm,
            "converged": False,
        }
        if row.sum() < min_total_count:
            records.append(blank)
            continue
        try:
            fit = bbreg(row, totals_array, formula, data, **kwargs)
        except Exception:
            # A guide that cannot be fitted is a data property, not a bug: the
            # screen continues and the guide is reported as unconverged, which
            # is how the R version behaves.
            records.append(blank)
            continue
        entry = fit.coefficient_table.loc[term]
        records.append(
            {
                "estimate": float(entry["estimate"]),
                "std_error": float(entry["std_error"]),
                "t_value": float(entry["t_value"]),
                "df": float(entry["df"]),
                "p_value": float(entry["p_value"]),
                "rho": fit.rho,
                "pearson_null": fit.pearson_null,
                "mean_cpm": mean_cpm,
                "converged": fit.converged,
            }
        )

    result = pd.DataFrame.from_records(records)
    result.insert(0, "guide", guide)
    if gene is not None:
        result.insert(0, "gene", list(gene))
    result["fdr"] = benjamini_hochberg(result["p_value"].to_numpy())
    return result


def bb_calibrate_controls(
    result: pd.DataFrame,
    control: Sequence[bool] | np.ndarray,
    alpha: float = 0.05,
    min_controls: int = 20,
    min_scale: float = 1.0,
    method: str = "tail_quantile",
) -> pd.DataFrame:
    """Rescale guide-level t tests against negative-control guides.

    Estimates a one-parameter empirical-null scale and divides every t
    statistic by it. ``min_scale`` defaults to 1 so that calibration can only
    make an already conservative analysis more conservative, never less.
    """
    if method not in {"tail_quantile", "qq_slope"}:
        raise ValueError("`method` must be 'tail_quantile' or 'qq_slope'.")
    if not (0 < alpha < 0.5):
        raise ValueError("`alpha` must be one finite number between 0 and 0.5.")
    if min_controls < 2:
        raise ValueError("`min_controls` must be at least two.")
    if min_scale <= 0:
        raise ValueError("`min_scale` must be positive.")

    control = np.asarray(control, dtype=bool)
    t_value = result["t_value"].to_numpy(dtype=float)
    df = result["df"].to_numpy(dtype=float)
    valid = control & np.isfinite(t_value) & np.isfinite(df)
    if int(valid.sum()) < int(min_controls):
        raise ValueError(f"At least {int(min_controls)} finite negative-control statistics are required.")

    control_df = np.unique(df[valid])
    if control_df.size != 1 or control_df[0] <= 0:
        raise ValueError("Finite negative-control guides must share one positive `df`.")
    control_df = float(control_df[0])

    absolute = np.abs(t_value[valid])
    if method == "tail_quantile":
        # One order statistic in the far tail. Simple, but high variance with
        # only a few hundred controls.
        empirical_cutoff = float(_quantile_type8(absolute, 1 - alpha))
        reference_cutoff = float(stats.t.ppf(1 - alpha / 2, df=control_df))
        ratio = empirical_cutoff / reference_cutoff
    else:
        # Slope of the control QQ plot against the t reference, through the
        # origin, over a band where many order statistics contribute.
        probabilities = np.round(np.arange(0.50, 0.9501, 0.01), 10)
        empirical = _quantile_type8(absolute, probabilities)
        reference = stats.t.ppf((1 + probabilities) / 2, df=control_df)
        ratio = float(np.sum(empirical * reference) / np.sum(reference**2))

    scale = max(min_scale, ratio)

    out = result.copy()
    out["raw_std_error"] = out["std_error"]
    out["raw_t_value"] = out["t_value"]
    out["raw_p_value"] = out["p_value"]
    out["raw_fdr"] = out["fdr"]
    out["std_error"] = out["std_error"] * scale
    out["t_value"] = out["t_value"] / scale
    out["p_value"] = 2 * stats.t.cdf(-np.abs(out["t_value"].to_numpy()), df=out["df"].to_numpy())
    out["fdr"] = benjamini_hochberg(out["p_value"].to_numpy())
    out.attrs["control_scale"] = scale
    out.attrs["control_alpha"] = alpha
    return out


def bb_gene_original(result: pd.DataFrame, min_guides: int = 1) -> pd.DataFrame:
    """Historical signed-z guide-to-gene aggregation.

    Combines each gene's guide statistics as a signed Stouffer z. This is the
    ``original`` statistic; the three later ones in the R implementation are
    not ported.
    """
    for column in ("gene", "estimate", "p_value"):
        if column not in result.columns:
            raise ValueError(
                "`result` must contain guide-level `gene`, `estimate`, and `p_value` columns."
            )
    if min_guides < 1:
        raise ValueError("`min_guides` must be one positive integer.")

    estimate = result["estimate"].to_numpy(dtype=float)
    p_value = result["p_value"].to_numpy(dtype=float)
    valid = np.isfinite(estimate) & np.isfinite(p_value) & (p_value >= 0) & (p_value <= 1)
    if "converged" in result.columns:
        converged = result["converged"].to_numpy()
        valid = valid & np.array([bool(c) if c is not None else False for c in converged])

    genes = result["gene"].to_numpy()
    rows = []
    # Sorted, to match R's `split()`, which orders groups by factor level.
    for gene_name in sorted(pd.unique(genes), key=str):
        all_index = np.flatnonzero(genes == gene_name)
        index = all_index[valid[all_index]]
        if index.size < min_guides:
            continue

        # Two-sided p to a signed z, floored at the smallest positive double so
        # that a p-value of exactly zero does not become an infinite z.
        signed_z = np.sign(estimate[index]) * stats.norm.isf(
            np.maximum(p_value[index] / 2, _DOUBLE_XMIN)
        )
        combined_z = float(np.sum(signed_z) / np.sqrt(index.size))
        gene_estimate = float(np.median(estimate[index]))
        rows.append(
            {
                "gene": gene_name,
                "n_guides": int(index.size),
                "estimate": gene_estimate,
                "statistic": combined_z,
                "guide_direction_agreement": (
                    np.nan
                    if gene_estimate == 0
                    else float(np.mean(np.sign(estimate[index]) == np.sign(gene_estimate)))
                ),
                "effect_statistic_sign_agreement": bool(
                    gene_estimate == 0 or np.sign(gene_estimate) == np.sign(combined_z)
                ),
                "p_value": float(2 * stats.norm.cdf(-abs(combined_z))),
                "converged_fraction": (
                    float(np.mean(result["converged"].to_numpy()[all_index].astype(float)))
                    if "converged" in result.columns
                    else np.nan
                ),
                "method": "original",
            }
        )

    if not rows:
        raise ValueError("No gene has enough finite guide results.")
    gene_result = pd.DataFrame(rows)
    gene_result["fdr"] = benjamini_hochberg(gene_result["p_value"].to_numpy())
    return gene_result
