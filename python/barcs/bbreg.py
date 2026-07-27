r"""Beta-binomial regression with t-based coefficient tests.

A thin rpy2 wrapper around ``R/bbreg.R``. Every function here calls straight
into that file; nothing in this module reimplements any of its numerics. It
fits the mean model

.. math::

    \operatorname{logit}(\mu_i) = x_i' \beta

under the beta-binomial variance

.. math::

    \operatorname{Var}(K_i) = n_i \mu_i (1 - \mu_i)
                              \{1 + (n_i - 1) \rho\}.

Why a wrapper rather than a port
---------------------------------
An earlier version of this module was a from-scratch NumPy/SciPy
reimplementation of ``R/bbreg.R``, kept in step with the original by a test
suite that required the two to agree to eight significant figures. That is one
more place for the two to drift, and it only covered five of the ten public
functions below -- the ones judged worth reimplementing. This module instead
sources ``R/bbreg.R`` once and calls its functions through rpy2, so there is
exactly one implementation and every exported R function is available here,
not just the ones a port happened to cover.

The cost is what the R implementation always cost: an R installation with
Rcpp/RcppArmadillo. Both are already repository dependencies declared in
``pixi.toml`` -- this module does not add a new requirement so much as stop
pretending the Python side doesn't have one.

Read ``R/bbreg.R`` for the statistical reasoning; it is the reference text and
is not repeated here. What lives in this module is the marshalling between
pandas/NumPy and R: data frames, formulas, matrices, and the handful of named
R attributes (such as ``bb_calibrate_controls``'s ``control_scale``) that get
attached to the returned DataFrame's ``.attrs``.

What changed for callers
-------------------------
``bb_screen``'s ``min_total_count`` now defaults to R's ``10`` rather than the
port's ``1`` -- the two were never supposed to disagree here, and now they
cannot. ``bbreg``/``bb_screen`` no longer accept a pre-built design matrix in
place of a formula string; R's formula machinery is now the only path, so
pass ``"~ dose + batch"`` rather than an already-expanded matrix.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter

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

# Sourced once, into R's global environment, the first time this module is
# imported. No CB2 package is loaded, so the C++ kernels in `R/bbreg.R` are
# never found and the pure-R fallback runs -- the same reference path
# `python/tests/reference_fit.R` exercises directly.
_R_SCRIPT = Path(__file__).resolve().parents[2] / "R" / "bbreg.R"
# Order matters, and pandas2ri must come last. Both converters claim R's
# `data.frame`: numpy2ri turns one into a structured `recarray`, pandas2ri into
# a DataFrame, and whichever is added later wins. Every function below returns
# a data frame, so pandas has to be the one that does.
_CONVERTER = robjects.default_converter + numpy2ri.converter + pandas2ri.converter
with localconverter(_CONVERTER):
    robjects.r["source"](str(_R_SCRIPT))
_R = robjects.globalenv


def _r_call(name: str, *args: Any, **kwargs: Any) -> Any:
    """Call one function sourced from ``R/bbreg.R``, converting both ways."""
    with localconverter(_CONVERTER):
        return _R[name](*args, **kwargs)


def _to_r(value: Any) -> Any:
    """Python -> R, passing anything that is already an R object straight through."""
    if isinstance(value, robjects.rinterface.Sexp):
        return value
    with localconverter(_CONVERTER):
        return robjects.conversion.get_conversion().py2rpy(value)


def _r_call_raw(name: str, *args: Any, **kwargs: Any) -> Any:
    """Call R, converting the arguments but returning the R object untouched.

    :func:`bbreg` needs this. Its result is an R list whose fields are pulled
    out one at a time by :meth:`BBRegResult._from_r`; if the return value were
    converted automatically it would arrive as a plain Python container with no
    ``.rx2`` to pull them out with.
    """
    return _R[name](
        *(_to_r(argument) for argument in args),
        **{key: _to_r(value) for key, value in kwargs.items()},
    )


def _to_py(robj: Any) -> Any:
    with localconverter(_CONVERTER):
        converted = robjects.conversion.get_conversion().rpy2py(robj)
    if isinstance(converted, pd.DataFrame):
        # R data frames carry row names, and for every frame BARCS returns
        # those are just the default "1", "2", ... Left alone they arrive as a
        # string index, so `screen["guide"][0]` raises KeyError and positional
        # code has to reach for `.iloc` everywhere. None of these frames uses
        # row names to carry information, so drop them for a RangeIndex.
        converted = converted.reset_index(drop=True)
    return converted


def _r_attr(robj: Any, name: str) -> Any:
    """Read one ``attr(robj, name)``, or ``None`` if R has not set it."""
    value = robjects.r["attr"](robj, name)
    if value is robjects.NULL:
        return None
    converted = _to_py(value)
    if hasattr(converted, "__len__") and not isinstance(converted, str) and len(converted) == 1:
        return converted[0]
    return converted


def _dataframe_with_attrs(robj: Any, attr_names: Sequence[str]) -> pd.DataFrame:
    """Convert an R data frame to pandas, carrying named R attributes into ``.attrs``.

    ``robj`` must be an unconverted R object -- attributes live on the R side
    and are lost the moment the frame becomes a DataFrame, so callers reach
    for :func:`_r_call_raw` rather than :func:`_r_call`.
    """
    frame = _to_py(robj)
    for name in attr_names:
        frame.attrs[name] = _r_attr(robj, name)
    return frame


def _matrix_to_frame(robj: Any) -> pd.DataFrame:
    """Convert an R matrix with row/column names to a labelled pandas DataFrame."""
    array = np.asarray(_to_py(robj))
    rownames = robjects.r["rownames"](robj)
    colnames = robjects.r["colnames"](robj)
    index = list(rownames) if rownames is not robjects.NULL else None
    columns = list(colnames) if colnames is not robjects.NULL else None
    return pd.DataFrame(array, index=index, columns=columns)


def _r_formula(formula: str | robjects.Formula) -> robjects.Formula:
    if isinstance(formula, robjects.Formula):
        return formula
    return robjects.Formula(str(formula))


def _named_r_vector(values: dict[str, float]) -> robjects.FloatVector:
    vector = robjects.FloatVector([float(v) for v in values.values()])
    vector.names = robjects.StrVector(list(values.keys()))
    return vector


def _optional_bool_vector(control: Sequence[bool] | np.ndarray | None) -> Any:
    if control is None:
        return robjects.NULL
    return robjects.BoolVector(np.asarray(control, dtype=bool))


def benjamini_hochberg(p_value: Sequence[float] | np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values, via R's ``p.adjust(method = "BH")``."""
    with localconverter(_CONVERTER):
        adjusted = robjects.r["p.adjust"](np.asarray(p_value, dtype=float), method="BH")
        return np.asarray(_to_py(adjusted))


@dataclass
class BBRegResult:
    """A fitted single-guide beta-binomial regression (R's ``bbreg`` object)."""

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
    # The underlying R `bbreg` object, kept so `bb_contrast` can pass this fit
    # straight back into R rather than refitting from the extracted fields.
    r_object: Any = field(repr=False, compare=False, default=None)

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

    @classmethod
    def _from_r(cls, r_fit: Any) -> "BBRegResult":
        table = _matrix_to_frame(r_fit.rx2("coefficient_table"))
        return cls(
            coefficients=pd.Series(_to_py(r_fit.rx2("coefficients")), index=table.index),
            coefficient_table=table,
            covariance=np.asarray(_to_py(r_fit.rx2("covariance"))),
            fitted_values=np.asarray(_to_py(r_fit.rx2("fitted.values"))),
            linear_predictors=np.asarray(_to_py(r_fit.rx2("linear.predictors"))),
            residuals=np.asarray(_to_py(r_fit.rx2("residuals"))),
            pearson=float(_to_py(r_fit.rx2("pearson"))[0]),
            pearson_null=float(_to_py(r_fit.rx2("pearson_null"))[0]),
            rho=float(_to_py(r_fit.rx2("rho"))[0]),
            scale=float(_to_py(r_fit.rx2("scale"))[0]),
            dispersion_boundary=bool(_to_py(r_fit.rx2("dispersion_boundary"))[0]),
            df_residual=int(_to_py(r_fit.rx2("df.residual"))[0]),
            rank=int(_to_py(r_fit.rx2("rank"))[0]),
            count=np.asarray(_to_py(r_fit.rx2("count"))),
            total=np.asarray(_to_py(r_fit.rx2("total"))),
            design=np.asarray(_to_py(r_fit.rx2("design"))),
            weights=np.asarray(_to_py(r_fit.rx2("weights"))),
            converged=bool(_to_py(r_fit.rx2("converged"))[0]),
            iterations=int(_to_py(r_fit.rx2("iterations"))[0]),
            r_object=r_fit,
        )


def bbreg(
    count: Sequence[float] | np.ndarray,
    total: Sequence[float] | np.ndarray,
    formula: str,
    data: pd.DataFrame,
    maxit: int = 100,
    tolerance: float = 1e-8,
    mu_bound: float = 1e-8,
) -> BBRegResult:
    """Fit beta-binomial regression for one guide. See ``R/bbreg.R::bbreg``."""
    r_fit = _r_call_raw(
        "bbreg",
        np.asarray(count, dtype=float),
        np.asarray(total, dtype=float),
        _r_formula(formula),
        data,
        maxit=maxit,
        tolerance=tolerance,
        mu_bound=mu_bound,
    )
    return BBRegResult._from_r(r_fit)


def bb_contrast(
    fit: BBRegResult,
    contrast: Sequence[float] | dict[str, float] | np.ndarray,
    null: float = 0.0,
) -> pd.DataFrame:
    """Test a linear contrast of the fitted coefficients. See ``R/bbreg.R::bb_contrast``."""
    r_contrast = (
        _named_r_vector(contrast)
        if isinstance(contrast, dict)
        else robjects.FloatVector(np.asarray(contrast, dtype=float))
    )
    r_result = _r_call("bb_contrast", fit.r_object, r_contrast, null=float(null))
    return _to_py(r_result)


def bb_screen(
    counts: np.ndarray | pd.DataFrame,
    data: pd.DataFrame,
    formula: str,
    term: str,
    totals: Sequence[float] | None = None,
    guide: Sequence[str] | None = None,
    gene: Sequence[str] | None = None,
    min_total_count: int = 10,
    ncores: int = 1,
    **kwargs: Any,
) -> pd.DataFrame:
    """Apply :func:`bbreg` guide by guide. See ``R/bbreg.R::bb_screen``.

    Returns one row per guide with the reported term's estimate, standard
    error, t statistic, degrees of freedom, p-value, BH-adjusted FDR,
    dispersion, and mean CPM.
    """
    if isinstance(counts, pd.DataFrame):
        if guide is None:
            guide = list(counts.index)
        counts = counts.to_numpy(dtype=float)
    counts_matrix = np.asarray(counts, dtype=float)

    call_kwargs: dict[str, Any] = {"min_total_count": min_total_count, "ncores": ncores}
    if totals is not None:
        call_kwargs["totals"] = np.asarray(totals, dtype=float)
    if guide is not None:
        call_kwargs["guide"] = robjects.StrVector([str(g) for g in guide])
    if gene is not None:
        call_kwargs["gene"] = robjects.StrVector([str(g) for g in gene])

    r_result = _r_call(
        "bb_screen", counts_matrix, data, _r_formula(formula), term, **call_kwargs, **kwargs
    )
    return _to_py(r_result)


def bb_calibrate_controls(
    result: pd.DataFrame,
    control: Sequence[bool] | np.ndarray,
    alpha: float = 0.05,
    min_controls: int = 20,
    min_scale: float = 1.0,
    method: str = "tail_quantile",
) -> pd.DataFrame:
    """Rescale guide-level t tests against negative-control guides.

    See ``R/bbreg.R::bb_calibrate_controls``. Returns ``result`` with
    recalibrated standard errors, t statistics, p-values, and FDR; the raw
    columns are kept with a ``raw_`` prefix, and the estimated scale and
    alpha are attached to the returned frame's ``.attrs``.
    """
    r_result = _r_call_raw(
        "bb_calibrate_controls",
        result,
        robjects.BoolVector(np.asarray(control, dtype=bool)),
        alpha=alpha,
        min_controls=min_controls,
        min_scale=min_scale,
        method=method,
    )
    return _dataframe_with_attrs(r_result, ["control_scale", "control_alpha"])


def bb_moderate_dispersion(
    result: pd.DataFrame,
    trend: bool = True,
    one_way: bool = False,
    borrow_df: bool = True,
    span: float = 0.5,
    min_guides: int = 50,
) -> pd.DataFrame:
    """Shrink guide-level dispersion toward a library-wide trend.

    See ``R/bbreg.R::bb_moderate_dispersion``. The estimated prior degrees of
    freedom and total reference degrees of freedom are attached to the
    returned frame's ``.attrs``.
    """
    r_result = _r_call_raw(
        "bb_moderate_dispersion",
        result,
        trend=bool(trend),
        one_way=bool(one_way),
        borrow_df=bool(borrow_df),
        span=span,
        min_guides=min_guides,
    )
    return _dataframe_with_attrs(r_result, ["prior_df", "df_total"])


def bb_gene_original(result: pd.DataFrame, min_guides: int = 1) -> pd.DataFrame:
    """The historical signed-z guide-to-gene statistic. See ``R/bbreg.R::bb_gene_original``."""
    r_result = _r_call("bb_gene_original", result, min_guides=min_guides)
    return _to_py(r_result)


def bb_gene_normal(
    result: pd.DataFrame,
    min_guides: int = 3,
    reference: str = "student_t",
) -> pd.DataFrame:
    """Exchangeable-normal guide-to-gene test. See ``R/bbreg.R::bb_gene_normal``."""
    r_result = _r_call_raw("bb_gene_normal", result, min_guides=min_guides, reference=reference)
    return _dataframe_with_attrs(r_result, ["reference", "null_assumption"])


def bb_gene_consistency(
    result: pd.DataFrame,
    control: Sequence[bool] | np.ndarray | None = None,
    min_guides: int = 3,
    alpha: float = 0.05,
    min_control_genes: int = 10,
    min_scale: float = 1.0,
) -> pd.DataFrame:
    """Empirical-null guide-consistency test. See ``R/bbreg.R::bb_gene_consistency``."""
    r_result = _r_call(
        "bb_gene_consistency",
        result,
        control=_optional_bool_vector(control),
        min_guides=min_guides,
        alpha=alpha,
        min_control_genes=min_control_genes,
        min_scale=min_scale,
    )
    return _dataframe_with_attrs(
        r_result,
        [
            "null_center",
            "null_scale",
            "global_scale",
            "control_scale",
            "control_genes",
            "null_assumption",
        ],
    )


def bb_gene_partial_pool(
    result: pd.DataFrame,
    control: Sequence[bool] | np.ndarray | None = None,
    min_guides: int = 2,
    alpha: float = 0.05,
    min_control_genes: int = 10,
    min_scale: float = 1.0,
) -> pd.DataFrame:
    """Random-effects guide pooling. See ``R/bbreg.R::bb_gene_partial_pool``."""
    r_result = _r_call_raw(
        "bb_gene_partial_pool",
        result,
        control=_optional_bool_vector(control),
        min_guides=min_guides,
        alpha=alpha,
        min_control_genes=min_control_genes,
        min_scale=min_scale,
    )
    return _dataframe_with_attrs(r_result, ["heterogeneity_estimator", "null_assumption"])


def bb_gene_eb_moderate(
    result: pd.DataFrame,
    control: Sequence[bool] | np.ndarray | None = None,
    min_guides: int = 2,
    prior_df: float = 4,
    alpha: float = 0.05,
    min_control_genes: int = 10,
    min_scale: float = 1.0,
) -> pd.DataFrame:
    """Empirical-Bayes moderated guide pooling. See ``R/bbreg.R::bb_gene_eb_moderate``."""
    r_result = _r_call_raw(
        "bb_gene_eb_moderate",
        result,
        control=_optional_bool_vector(control),
        min_guides=min_guides,
        prior_df=prior_df,
        alpha=alpha,
        min_control_genes=min_control_genes,
        min_scale=min_scale,
    )
    return _dataframe_with_attrs(
        r_result, ["prior_tau2", "prior_df", "heterogeneity_estimator", "null_assumption"]
    )
