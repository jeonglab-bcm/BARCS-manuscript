# Locating the external MAGeCK installation.
#
# MAGeCK is deliberately not vendored into this repository. It is the external
# negative-binomial benchmark BARCS is measured against, and embedding a
# partial reimplementation would make that comparison worthless. So the
# benchmarks shell out to the official executable, and something has to decide
# where that executable is.
#
# Historically each script hardcoded `.venv/bin/mageck`, which was where an
# ad-hoc virtualenv happened to put it on one machine. That is now the last
# resort rather than the only option, because `pixi.toml` declares MAGeCK as a
# conda dependency and puts it on PATH.
#
# Resolution order, most specific first:
#
#   1. $MAGECK              -- an explicit override, for testing another build
#   2. mageck on PATH       -- what `pixi run` and `pixi shell` provide
#   3. .venv/bin/mageck     -- the legacy ad-hoc virtualenv
#
# Source this file, then call `mageck_executable()`.

# Every committed MAGeCK number was produced with this version. Its gene
# p-value is a permutation tail probability, so a different release can shift
# results without anything in this repository changing. A mismatch warns
# rather than stops: reproducing an old number needs this version, but trying
# a new one is a legitimate thing to want to do.
MAGECK_PINNED_VERSION <- "0.5.9.5"

.mageck_legacy_venv <- file.path(".venv", "bin", "mageck")

#' Path to the official MAGeCK executable
#'
#' @param required Stop with an actionable message when MAGeCK cannot be
#'   found. Set FALSE to test for availability instead.
#' @return Path to the executable, or "" when absent and `required` is FALSE.
mageck_executable <- function(required = TRUE) {
  override <- Sys.getenv("MAGECK", unset = "")
  if (nzchar(override)) {
    if (!file.exists(override)) {
      stop("$MAGECK is set to '", override, "', which does not exist.",
           call. = FALSE)
    }
    return(override)
  }

  on_path <- Sys.which("mageck")
  if (nzchar(on_path)) {
    return(unname(on_path))
  }

  if (file.exists(.mageck_legacy_venv)) {
    return(.mageck_legacy_venv)
  }

  if (!required) {
    return("")
  }
  stop(
    "Official MAGeCK ", MAGECK_PINNED_VERSION, " is required but was not ",
    "found.\n",
    "  The declared toolchain provides it:  pixi run doctor\n",
    "  Or point at an existing install:     MAGECK=/path/to/mageck Rscript ...",
    call. = FALSE
  )
}

#' The Python interpreter that can import MAGeCK's modules
#'
#' Two benchmarks do not invoke MAGeCK as a command. They run
#' `scripts/mageck_compat.py`, which imports MAGeCK's own entry point in order
#' to patch two NumPy behaviours it still depends on, and
#' `scripts/mageck_cnv_correct.py`, which calls its CNV normalizer directly.
#' Both must run under an interpreter that can see the MAGeCK installation,
#' which is the one alongside the executable -- not necessarily the first
#' `python3` on PATH.
mageck_python <- function(required = TRUE) {
  override <- Sys.getenv("MAGECK_PYTHON", unset = "")
  if (nzchar(override)) {
    return(override)
  }

  executable <- mageck_executable(required = required)
  if (!nzchar(executable)) {
    return("")
  }

  # MAGeCK is installed as a console script, so its interpreter is its
  # neighbour in the same bin directory under both pixi and a virtualenv.
  for (candidate in c("python3", "python")) {
    beside <- file.path(dirname(executable), candidate)
    if (file.exists(beside)) {
      return(beside)
    }
  }

  fallback <- Sys.which("python3")
  if (nzchar(fallback)) {
    return(unname(fallback))
  }
  if (!required) {
    return("")
  }
  stop("No Python interpreter found alongside MAGeCK at '", executable, "'.",
       call. = FALSE)
}

#' Installed MAGeCK version, as a string
mageck_version <- function() {
  raw <- tryCatch(
    system2(mageck_executable(), "--version", stdout = TRUE, stderr = TRUE),
    error = function(e) NA_character_
  )
  trimws(paste(raw, collapse = ""))
}

#' Warn when the installed MAGeCK is not the version the results were built on
mageck_check_version <- function() {
  found <- mageck_version()
  if (!is.na(found) && !identical(found, MAGECK_PINNED_VERSION)) {
    warning(
      "MAGeCK ", found, " is installed, but the committed results were ",
      "produced with ", MAGECK_PINNED_VERSION, ". Its permutation null makes ",
      "gene p-values version sensitive.",
      call. = FALSE
    )
  }
  invisible(found)
}

#' A PATH= string that puts MAGeCK first
#'
#' `system2(env = )` replaces rather than extends the environment, so callers
#' that need MAGeCK's own bin directory visible to a subprocess pass this.
mageck_path_env <- function() {
  paste0(
    "PATH=", normalizePath(dirname(mageck_executable())), ":",
    Sys.getenv("PATH")
  )
}
