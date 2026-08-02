# Make the BARCS package available to a manuscript script.
#
# BARCS used to live in this repository as a loose `R/bbreg.R` that every
# script sourced. It is now its own package, tracked here as the `BARCS/`
# submodule, so there is exactly one copy of the implementation and each
# manuscript revision records the package commit it was produced with.
#
# Every analysis script begins with:
#
#   source(file.path("R", "load_barcs.R"))
#
# from the repository root. That prefers an installed BARCS and otherwise
# loads the pinned submodule in place, so a clean checkout works with no
# install step.

local({
  minimum_version <- "0.1.0"

  loaded_from <- if (requireNamespace("BARCS", quietly = TRUE)) {
    suppressPackageStartupMessages(
      library("BARCS", character.only = TRUE)
    )
    paste0("installed package ", utils::packageVersion("BARCS"))
  } else if (file.exists(file.path("BARCS", "DESCRIPTION"))) {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop(
        "BARCS is not installed and `pkgload` is unavailable to load the ",
        "submodule. Install either one:\n",
        "  Rscript -e 'install.packages(\"pkgload\")'\n",
        "  R CMD INSTALL BARCS",
        call. = FALSE
      )
    }
    # `load_all()` compiles the RcppArmadillo kernels on first use, so the
    # submodule path is as fast as an installed package after that.
    pkgload::load_all("BARCS", quiet = TRUE, export_all = FALSE)
    paste0("BARCS/ submodule ", utils::packageVersion("BARCS"))
  } else {
    stop(
      "Cannot find BARCS. Run this script from the repository root with the ",
      "submodule checked out:\n",
      "  git submodule update --init --recursive",
      call. = FALSE
    )
  }

  if (utils::packageVersion("BARCS") < minimum_version) {
    stop(
      sprintf(
        "BARCS >= %s is required; found %s. Update the submodule with:\n%s",
        minimum_version, utils::packageVersion("BARCS"),
        "  git submodule update --remote BARCS"
      ),
      call. = FALSE
    )
  }

  message("BARCS loaded from ", loaded_from)
})
