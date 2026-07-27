#!/usr/bin/env Rscript

# Install the two R packages that exist only on GitHub.
#
# Everything else the analysis needs is a conda package declared in
# `pixi.toml`, which means `pixi.lock` pins its exact build. These two are the
# exception: neither is on CRAN, Bioconductor, or any conda channel, so they
# are pinned here instead, by commit.
#
#   CRISPhieRmix  A hierarchical mixture over the guides within a gene. It is
#                 the gene caller BARCS is compared against at genome scale,
#                 chosen over the RNA-seq count models because it was actually
#                 designed for pooled screens.
#
#   simCRISPR     The factorial screen simulator behind the interaction
#                 benchmark. Unlike CRISPulator it gives each guide a knockout
#                 effect, a treatment effect, and an interaction, which is what
#                 makes it a test of the coefficient BARCS exists to estimate.
#
# Pinning to a commit rather than to a branch matters for the same reason the
# Julia version is pinned: both packages feed simulated data or gene calls into
# committed metrics, and an upstream change would silently move those numbers.
#
# Packages install into the active pixi environment's own R library, so they
# disappear with `pixi clean` and never leak into a user-level library.
#
#     pixi run setup-r-github

github_packages <- list(
  list(
    name = "CRISPhieRmix",
    repo = "timydaley/CRISPhieRmix",
    # HEAD of master as of 2026-07-27; the repository was last touched in 2019.
    ref  = "e400f21e00bfdb2ec5e88fb6e1dc8b6e9fba1cb5"
  ),
  list(
    name = "simCRISPR",
    repo = "bachergroup/simCRISPR",
    # HEAD of main as of 2026-07-27.
    ref  = "677e95c5b463911d094a685398ccc7e610a354f2"
  )
)

# Install into the first writable library on the path, which inside a pixi
# task is the environment's own. Bail out rather than fall back to a
# user-level library, since that would put the package somewhere `pixi clean`
# cannot reach and make the environment non-reproducible.
target_lib <- .libPaths()[1]
if (file.access(target_lib, mode = 2) != 0) {
  stop(
    "R library '", target_lib, "' is not writable. Run this through pixi:\n",
    "    pixi run setup-r-github",
    call. = FALSE
  )
}

for (pkg in github_packages) {
  if (requireNamespace(pkg$name, quietly = TRUE)) {
    message("already installed: ", pkg$name)
    next
  }
  message("installing ", pkg$name, " from ", pkg$repo, "@", substr(pkg$ref, 1, 7))
  remotes::install_github(
    paste0(pkg$repo, "@", pkg$ref),
    lib = target_lib,
    # Dependencies are declared in pixi.toml so the lock file covers them;
    # letting remotes resolve them here would install unpinned copies on top.
    dependencies = FALSE,
    upgrade = "never",
    quiet = FALSE
  )
}

# Fail loudly now rather than three hundred lines into a benchmark.
missing <- Filter(
  function(p) !requireNamespace(p$name, quietly = TRUE),
  github_packages
)
if (length(missing) > 0) {
  stop(
    "failed to install: ",
    paste(vapply(missing, function(p) p$name, character(1)), collapse = ", "),
    call. = FALSE
  )
}

message("GitHub R packages ready in ", target_lib)
