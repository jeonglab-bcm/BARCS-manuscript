# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "formulaic>=1.2.2",
#     "marimo>=0.23.3",
#     "matplotlib>=3.11.1",
#     "numpy>=2.5.1",
#     "pandas>=3.0.5",
#     "scipy>=1.18.0",
# ]
# ///
"""Interactive walkthrough of BARCS.

A marimo notebook. Marimo stores notebooks as plain Python, which is why this
is a `.py` file and not a `.ipynb`: it diffs, it reviews, and it imports like
any other module.

    pixi run notebook        # edit
    pixi run notebook-run    # read-only

The notebook is reactive rather than sequential. Changing a slider re-runs
every cell that depends on it and nothing else, so there is no stale-state
problem and no need to remember which cells were run in what order -- which
matters here, because the point of the notebook is to move the simulation
parameters around and watch the inference respond.

What it is for: showing what BARCS does and where it breaks, on data whose
truth is known. It is not a benchmark. The benchmarks are the R scripts in
`examples/`, and `docs/repository-map.md` says which one answers what.
"""

import marimo

__generated_with = "0.23.4"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # BARCS: beta-binomial regression for CRISPR screens

    A pooled screen measures each guide as a count out of a library total.
    The obvious model for a count out of a total is the binomial, and it is
    wrong in a way that matters: it assumes the only variability is
    sequencing depth. Real replicate libraries differ by more than that, so
    the binomial understates the variance, and a test built on it calls
    noise significant.

    BARCS keeps the binomial sampling step and adds a beta-distributed
    guide proportion on top of it:

    $$\operatorname{Var}(K_i) = n_i \mu_i (1 - \mu_i)
      \bigl\{1 + (n_i - 1)\rho\bigr\}$$

    with $\rho$ estimated per guide. The mean model is an ordinary design
    matrix, $\operatorname{logit}(\mu_i) = x_i'\beta$, which is what lets
    a screen be analysed against a dose, a time course, or an adjusted
    contrast rather than only two groups.

    Everything below is simulated, so the true effects are known and the
    answers can be checked rather than admired.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## 1. Simulate a screen
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    # Defaults are a small but realistic dose-response screen. `rho` is the
    # control that matters most: at 0 the binomial would be correct and BARCS
    # has nothing to fix, and the naive comparison below collapses to a tie.
    n_genes = mo.ui.slider(20, 300, value=120, step=20, label="genes")
    guides_per_gene = mo.ui.slider(2, 8, value=4, label="guides per gene")
    active_fraction = mo.ui.slider(
        0.0, 0.5, value=0.2, step=0.05, label="fraction of genes with an effect"
    )
    effect_size = mo.ui.slider(
        0.0, 2.0, value=0.8, step=0.1, label="effect size (logit per dose unit)"
    )
    rho = mo.ui.slider(
        0.0, 0.01, value=0.002, step=0.0005, label="true rho (overdispersion)"
    )
    depth = mo.ui.slider(
        20_000, 200_000, value=80_000, step=20_000, label="median library size"
    )
    seed = mo.ui.number(1, 10_000, value=42, label="seed")

    mo.vstack(
        [
            mo.hstack([n_genes, guides_per_gene], justify="start"),
            mo.hstack([active_fraction, effect_size], justify="start"),
            mo.hstack([rho, depth], justify="start"),
            seed,
        ]
    )
    return (
        active_fraction,
        depth,
        effect_size,
        guides_per_gene,
        n_genes,
        rho,
        seed,
    )


@app.cell
def _(
    active_fraction,
    depth,
    effect_size,
    guides_per_gene,
    n_genes,
    np,
    pd,
    rho,
    seed,
):
    def simulate_screen(
        n_genes, guides_per_gene, active_fraction, effect_size, rho, depth, seed
    ):
        """Draw a dose-response screen from the model BARCS assumes.

        Five dose levels in three batches. Counts come from a beta-binomial:
        each guide has a latent proportion drawn from a Beta with the given
        intraclass correlation, and the observed count is a binomial draw from
        that. Drawing from the assumed model is the point -- it makes the truth
        exactly known, so the estimator can be checked rather than compared.
        """
        rng = np.random.default_rng(int(seed))

        dose = np.tile([-1.5, -0.75, 0.0, 0.75, 1.5], 3)
        batch = np.repeat(["b1", "b2", "b3"], 5)
        samples = pd.DataFrame({"dose": dose, "batch": pd.Categorical(batch)})

        totals = rng.integers(int(depth * 0.6), int(depth * 1.4), size=len(samples))
        totals = totals.astype(float)

        n_active = int(round(n_genes * active_fraction))
        gene_effect = np.zeros(n_genes)
        if n_active:
            signs = rng.choice([-1.0, 1.0], size=n_active)
            gene_effect[:n_active] = signs * effect_size

        gene = np.repeat([f"gene{i:04d}" for i in range(n_genes)], guides_per_gene)
        guide = np.array(
            [f"{name}_sg{i % guides_per_gene}" for i, name in enumerate(gene)]
        )
        guide_effect = np.repeat(gene_effect, guides_per_gene)
        n_guides = n_genes * guides_per_gene

        baseline = rng.uniform(-7.8, -6.5, size=n_guides)
        offsets = np.array([{"b1": 0.0, "b2": 0.18, "b3": -0.12}[b] for b in batch])

        counts = np.zeros((n_guides, len(samples)))
        for g in range(n_guides):
            eta = baseline[g] + guide_effect[g] * dose + offsets
            mu = 1.0 / (1.0 + np.exp(-eta))
            if rho > 0:
                precision = 1.0 / rho - 1.0
                proportion = rng.beta(mu * precision, (1 - mu) * precision)
            else:
                # rho = 0 is the binomial special case; the Beta degenerates.
                proportion = mu
            counts[g] = rng.binomial(totals.astype(int), proportion)

        truth = pd.DataFrame(
            {
                "gene": [f"gene{i:04d}" for i in range(n_genes)],
                "true_effect": gene_effect,
                "active": gene_effect != 0,
            }
        )
        return counts, samples, totals, gene, guide, truth

    counts, samples, totals, gene, guide, truth = simulate_screen(
        n_genes.value,
        guides_per_gene.value,
        active_fraction.value,
        effect_size.value,
        rho.value,
        depth.value,
        seed.value,
    )
    return counts, gene, guide, samples, totals, truth


@app.cell(hide_code=True)
def _(counts, mo, truth):
    mo.md(
        f"""
    Simulated **{counts.shape[0]:,} guides** over **{counts.shape[1]} samples**
    ({len(truth)} genes, {int(truth["active"].sum())} of them active).
    Median count per guide per sample: **{counts.mean():,.0f}**.

    The design is a continuous dose crossed with three batches:
    `~ dose + batch`. That is the whole reason for the regression --
    the dose column is a number, not a label, and the batch adjustment
    happens inside the same fit rather than as a preprocessing step.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## 2. One guide, in detail

    Before the screen, a single fit. This is what `bbreg()` returns for one
    guide: a coefficient table indexed by the design-matrix columns, plus
    the estimated `rho`.
    """)
    return


@app.cell(hide_code=True)
def _(guide, mo):
    guide_index = mo.ui.dropdown(
        options={name: i for i, name in enumerate(guide)},
        value=guide[0],
        label="guide",
        searchable=True,
        full_width=True,
    )
    guide_index
    return (guide_index,)


@app.cell
def _(
    bbreg,
    counts,
    guide,
    guide_index,
    guides_per_gene,
    mo,
    samples,
    totals,
    truth,
):
    # R does the fitting; `BBRegResult` is the R object's fields unpacked into
    # a dataclass, with the untouched R object kept on `.r_object` so
    # `bb_contrast` can pass it straight back rather than refitting.
    single_fit = bbreg(
        counts[guide_index.value],
        totals,
        "~ dose + batch",
        samples,
    )
    true_for_guide = truth["true_effect"][guide_index.value // guides_per_gene.value]

    mo.vstack(
        [
            mo.md(
                f"**{guide[guide_index.value]}** &nbsp;·&nbsp; "
                f"true dose effect: **{true_for_guide:+.2f}** &nbsp;·&nbsp; "
                f"estimated rho: `{single_fit.rho:.2e}` &nbsp;·&nbsp; "
                f"converged: `{single_fit.converged}`"
            ),
            single_fit.coefficient_table.round(6),
        ]
    )
    return (single_fit,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ### A contrast

    Coefficients answer "per one unit of dose". Questions that span more
    than one unit, or that combine coefficients, are contrasts.
    `bb_contrast()` takes the linear combination and carries the covariance
    through, which is not the same as scaling the coefficient's standard
    error by hand once covariates are correlated.
    """)
    return


@app.cell
def _(bb_contrast, single_fit):
    # Two dose steps rather than one.
    bb_contrast(single_fit, {"dose": 2.0}).round(6)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## 3. The whole screen, and why the binomial is not enough

    `bb_screen()` fits every guide and reports one coefficient. Alongside
    it, the same screen analysed with a plain binomial GLM -- the model that
    assumes sequencing depth is the only source of variation.

    Watch the **null rejection rate**. Under a correct model it should sit
    near 0.05. Drag `rho` toward zero and the two converge; push it up and
    the binomial's error rate climbs while BARCS stays put. That gap is the
    entire argument for the method.
    """)
    return


@app.cell
def _(bb_screen, counts, gene, guide, samples, sm_binomial, totals):
    screen = bb_screen(
        counts,
        samples,
        "~ dose + batch",
        term="dose",
        totals=totals,
        guide=list(guide),
        gene=list(gene),
    )
    naive = sm_binomial(counts, samples, totals)
    screen_with_naive = screen.assign(
        naive_p_value=naive["p_value"], naive_estimate=naive["estimate"]
    )
    return (screen_with_naive,)


@app.cell
def _(expit, norm, np, pd):
    def sm_binomial(counts, samples, totals):
        """Binomial IRLS per guide -- the deliberately misspecified comparison.

        Written out rather than pulled from statsmodels so the notebook has one
        fewer dependency, and so it is visibly the *same* design matrix and the
        *same* iteration as BARCS with the overdispersion term removed. The
        difference in the results below is then attributable to that term
        alone.
        """
        from formulaic import model_matrix

        x = np.asarray(model_matrix("~ dose + batch", samples), dtype=float)
        dose_column = [
            i
            for i, name in enumerate(model_matrix("~ dose + batch", samples).columns)
            if name == "dose"
        ][0]
        df_residual = x.shape[0] - x.shape[1]

        estimates = np.full(counts.shape[0], np.nan)
        p_values = np.full(counts.shape[0], np.nan)
        for g in range(counts.shape[0]):
            proportion = counts[g] / totals
            mu = (counts[g] + 0.5) / (totals + 1.0)
            beta = np.zeros(x.shape[1])
            for _ in range(50):
                eta = np.log(mu / (1 - mu))
                weight = totals * mu * (1 - mu)
                working = eta + (proportion - mu) / (mu * (1 - mu))
                try:
                    beta_new = np.linalg.solve(
                        x.T @ (weight[:, None] * x), x.T @ (weight * working)
                    )
                except np.linalg.LinAlgError:
                    beta_new = beta
                    break
                if np.max(np.abs(beta_new - beta)) < 1e-10:
                    beta = beta_new
                    break
                beta = beta_new
                mu = np.clip(expit(x @ beta), 1e-10, 1 - 1e-10)

            try:
                weight = totals * mu * (1 - mu)
                covariance = np.linalg.inv(x.T @ (weight[:, None] * x))
            except np.linalg.LinAlgError:
                continue
            standard_error = np.sqrt(covariance[dose_column, dose_column])
            z = beta[dose_column] / standard_error
            estimates[g] = beta[dose_column]
            # A z reference, not a t: this is exactly the over-confidence being
            # illustrated, so it is not quietly corrected here.
            p_values[g] = 2 * norm.sf(abs(z))

        return pd.DataFrame({"estimate": estimates, "p_value": p_values})

    return (sm_binomial,)


@app.cell
def _(bb_gene_original, mo, pd, screen_with_naive, truth):
    genes = bb_gene_original(screen_with_naive).merge(truth, on="gene", how="left")

    null_genes = genes[~genes["active"]]
    active_genes = genes[genes["active"]]

    # The naive comparison is made at the guide level, where the binomial's
    # over-confidence originates, before any gene-level pooling smooths it.
    guide_truth = screen_with_naive.merge(truth, on="gene", how="left")
    guide_null = guide_truth[~guide_truth["active"]]

    summary = pd.DataFrame(
        {
            "model": ["BARCS (beta-binomial t)", "naive binomial z"],
            "null rejection rate at p<0.05": [
                float((guide_null["p_value"] < 0.05).mean()),
                float((guide_null["naive_p_value"] < 0.05).mean()),
            ],
            "guides tested": [len(guide_null), len(guide_null)],
        }
    )

    mo.vstack(
        [
            mo.md("### Guide-level type I error on genes with no effect"),
            mo.md(
                "A correctly calibrated test rejects 5% of true nulls at "
                "p<0.05. At the default settings the naive binomial rejects "
                "roughly **85%** of them -- it is not a slightly optimistic "
                "test, it is uninformative.\n\n"
                "BARCS lands near 9%, which is much better and still not 5%. "
                "That residual is a real limitation rather than a rounding "
                "error: the covariance treats each guide's dispersion as a "
                "fixed plug-in value, so the t reference does not carry the "
                "uncertainty in having *estimated* rho from eleven residual "
                "degrees of freedom. Section 5 is what a real screen does "
                "about it. Set `rho` to 0 and the two models converge, which "
                "is the check that the difference comes from overdispersion "
                "and not from something else in the comparison."
            ),
            summary.round(4),
        ]
    )
    return genes, guide_null, guide_truth


@app.cell(hide_code=True)
def _(genes, mo):
    called = genes[genes["fdr"] < 0.10]
    true_positive = int((called["active"]).sum())
    false_positive = int((~called["active"]).sum())
    realized = false_positive / max(len(called), 1)
    recall = true_positive / max(int(genes["active"].sum()), 1)

    mo.md(
        f"""
        ### Gene-level calls at FDR 0.10

        | | |
        |---|---|
        | genes called | **{len(called)}** |
        | truly active among them | **{true_positive}** |
        | realized false-discovery proportion | **{realized:.3f}** (requested 0.10) |
        | recall | **{recall:.3f}** |

        Realized FDP is the number to read, not recall. A method that calls
        everything has perfect recall; the question is whether the error rate
        it delivers is the one it advertised. Raise `rho` and watch what
        happens to the realized proportion.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## 4. What the fit looks like

    Left: estimated against true gene effect. Right: the p-value
    distribution for genes with no effect, which should be flat -- a spike
    near zero is exactly the miscalibration the table above counts.
    """)
    return


@app.cell(hide_code=True)
def _(genes, guide_null, np, plt):
    figure, (left, right) = plt.subplots(1, 2, figsize=(11, 4.2))

    left.axhline(0, color="#BBBBBB", lw=0.8)
    left.axvline(0, color="#BBBBBB", lw=0.8)
    left.scatter(
        genes["true_effect"],
        genes["estimate"],
        s=18,
        alpha=0.6,
        color="#0072B2",
        edgecolor="none",
    )
    limit = float(np.nanmax(np.abs(genes["true_effect"]))) * 1.3 + 0.1
    left.plot([-limit, limit], [-limit, limit], color="#D55E00", lw=1.2, ls="--")
    left.set_xlabel("true gene effect")
    left.set_ylabel("estimated gene effect")
    left.set_title("Effect recovery (dashed = identity)")

    right.hist(
        guide_null["p_value"].dropna(),
        bins=20,
        range=(0, 1),
        color="#0072B2",
        alpha=0.75,
        label="BARCS",
    )
    right.hist(
        guide_null["naive_p_value"].dropna(),
        bins=20,
        range=(0, 1),
        histtype="step",
        color="#D55E00",
        lw=1.6,
        label="naive binomial",
    )
    right.axhline(
        len(guide_null.dropna(subset=["p_value"])) / 20,
        color="#555555",
        ls=":",
        lw=1.2,
        label="uniform",
    )
    right.set_xlabel("p-value on guides with no true effect")
    right.set_ylabel("guides")
    right.set_title("Null calibration")
    right.legend(frameon=False, fontsize=9)

    figure.tight_layout()
    figure
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## 5. Negative-control calibration

    Real screens carry non-targeting guides whose true effect is zero. If
    their statistics are wider than the reference distribution says they
    should be, the whole screen is over-confident and can be rescaled to
    match — a one-parameter empirical null.

    Here the true nulls stand in for those controls. `min_scale` defaults
    to 1, so calibration can only ever make the analysis more conservative,
    never less.
    """)
    return


@app.cell
def _(bb_calibrate_controls, guide_truth, mo, pd, screen_with_naive):
    control = (~guide_truth["active"]).to_numpy()

    rows = []
    calibrated_frames = {}
    for method in ("tail_quantile", "qq_slope"):
        calibrated = bb_calibrate_controls(
            screen_with_naive, control=control, method=method, min_controls=10
        )
        calibrated_frames[method] = calibrated
        rows.append(
            {
                "method": method,
                "estimated null scale": calibrated.attrs["control_scale"],
                "null rejections at p<0.05": float(
                    (calibrated.loc[control, "p_value"] < 0.05).mean()
                ),
            }
        )

    mo.vstack(
        [
            mo.md(
                "A scale near 1 means the model was already calibrated and "
                "there is nothing to correct. `qq_slope` uses the whole "
                "0.50–0.95 band rather than a single order statistic, so it is "
                "the steadier of the two."
            ),
            pd.DataFrame(rows).round(4),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## 6. The four guide-to-gene statistics

    Guides have to be combined into genes somehow, and the choice is not
    incidental -- it is what the manuscript's later comparisons are about.
    All four are reachable here because `barcs` calls R rather than
    reimplementing it.

    | | |
    |---|---|
    | `original` | signed-z aggregation, the historical default |
    | `normal` | exchangeable normal guide coefficients |
    | `partial_pool` | random-effects partial pooling |
    | `eb_moderate` | empirical-Bayes moderation of guide heterogeneity |

    They receive identical guide-level fits, so any difference below is
    attributable to the pooling rule and nothing else.
    """)
    return


@app.cell
def _(
    bb_gene_eb_moderate,
    bb_gene_normal,
    bb_gene_original,
    bb_gene_partial_pool,
    mo,
    pd,
    screen_with_naive,
    truth,
):
    _statistics = {
        "original": bb_gene_original,
        "normal": bb_gene_normal,
        "partial_pool": bb_gene_partial_pool,
        "eb_moderate": bb_gene_eb_moderate,
    }

    _rows = []
    for _label, _function in _statistics.items():
        try:
            _scored = _function(screen_with_naive).merge(truth, on="gene", how="left")
        except Exception as _error:
            # Some statistics need a minimum number of guides per gene, so a
            # small slider setting can legitimately put one out of range. Say
            # which, rather than blanking the whole cell.
            _rows.append({"statistic": _label, "note": str(_error)[:70]})
            continue
        _called = _scored[_scored["fdr"] < 0.10]
        _rows.append(
            {
                "statistic": _label,
                "genes called": len(_called),
                "true positives": int(_called["active"].sum()),
                "realized FDP": (~_called["active"]).sum() / max(len(_called), 1),
                "recall": int(_called["active"].sum())
                / max(int(_scored["active"].sum()), 1),
            }
        )

    mo.vstack(
        [
            mo.md("### At gene FDR 0.10"),
            pd.DataFrame(_rows).round(4),
            mo.md(
                "Read realized FDP against the requested 0.10 first, then "
                "recall. A statistic that calls more genes has not necessarily "
                "done better -- it may simply be delivering a higher error rate "
                "than it advertised."
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## 7. Where this stops

    This notebook and the `barcs` Python package cover the guide-level
    inference path: `bbreg`, `bb_contrast`, `bb_screen`,
    `bb_calibrate_controls`, and all four guide-to-gene statistics.

    Everything above ran in R. `barcs` is an rpy2 wrapper around
    `R/bbreg.R` and reimplements none of it, so nothing here can drift away
    from what the benchmarks and the manuscript used, and the whole API is
    reachable rather than the subset somebody got round to porting.

    `pixi run test-python` checks the marshalling, including a run of this
    same pipeline against a native `Rscript` over one fixture, required to
    agree to 1e-12.

    This notebook is a demonstration, not a benchmark. The benchmarks are the
    R scripts in `examples/`, and `docs/repository-map.md` says which one
    answers what.
    """)
    return


@app.cell
def _():
    import sys
    from pathlib import Path

    # The package lives beside this notebook rather than being pip-installed,
    # so that editing `python/barcs/` is picked up on the next reactive run.
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from scipy.special import expit
    from scipy.stats import norm

    # `barcs` is an rpy2 wrapper around `R/bbreg.R`, so everything below runs
    # the same code the benchmarks and the manuscript used -- including the
    # three later gene statistics, which exist only in R.
    from barcs import (
        bb_calibrate_controls,
        bb_contrast,
        bb_gene_eb_moderate,
        bb_gene_normal,
        bb_gene_original,
        bb_gene_partial_pool,
        bb_screen,
        bbreg,
    )

    return (
        bb_calibrate_controls,
        bb_contrast,
        bb_gene_eb_moderate,
        bb_gene_normal,
        bb_gene_original,
        bb_gene_partial_pool,
        bb_screen,
        bbreg,
        expit,
        mo,
        norm,
        np,
        pd,
        plt,
    )


if __name__ == "__main__":
    app.run()
