#!/usr/bin/env bash
#
# Report whether this machine can actually run each part of the analysis.
#
#     pixi run doctor
#
# The repository drives four separate toolchains -- R for the method and the
# benchmarks, Julia for the CRISPulator simulations, Python for MAGeCK, and TeX
# for the manuscript -- and a missing piece of any one of them tends to surface
# hundreds of lines into a long-running script. This checks all of them up
# front and prints one line per capability.
#
# Exit status is 0 when everything required is present, 1 otherwise. Checks
# that are only needed for part of the work are reported as WARN and do not
# fail the run.

set -uo pipefail  # deliberately not -e: a failing check must be reported, not fatal

ROOT="${PIXI_PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "${ROOT}"

pass_count=0
warn_count=0
fail_count=0

pass() { printf '  \033[32mok\033[0m    %-22s %s\n' "$1" "${2:-}"; pass_count=$((pass_count + 1)); }
warn() { printf '  \033[33mwarn\033[0m  %-22s %s\n' "$1" "${2:-}"; warn_count=$((warn_count + 1)); }
fail() { printf '  \033[31mFAIL\033[0m  %-22s %s\n' "$1" "${2:-}"; fail_count=$((fail_count + 1)); }
section() { printf '\n\033[1m%s\033[0m\n' "$1"; }

# ---------------------------------------------------------------------------
section "R"
# ---------------------------------------------------------------------------

if ! command -v Rscript >/dev/null 2>&1; then
  fail "Rscript" "not on PATH -- run inside 'pixi shell' or via 'pixi run'"
else
  pass "Rscript" "$(Rscript -e 'cat(R.version.string)' 2>/dev/null)"

  # One R process checks every package, because R startup dominates the cost.
  # Each package is reported with what it is needed for, so a missing one says
  # which benchmarks are blocked rather than just that something is absent.
  #
  # The results go to a temp file instead of down a pipe: a `while read` loop
  # on the right of a pipe runs in a subshell, and the pass/fail counters it
  # incremented would be discarded when that subshell exited.
  r_report="$(mktemp)"
  Rscript --vanilla - >"${r_report}" 2>/dev/null <<'EOF'
packages <- list(
  c("CB2",          "the package under test (submodule)"),
  c("Rcpp",         "CB2 compiled core"),
  c("DESeq2",       "count-model comparator, CRISPhieRmix input"),
  c("edgeR",        "count-model comparator"),
  c("limma",        "count-model comparator"),
  c("CRISPhieRmix", "genome-scale gene-caller comparator"),
  c("simCRISPR",    "factorial interaction benchmark"),
  c("readxl",       "published supplementary tables"),
  c("rhdf5",        "Chronos HDF5 releases"),
  c("testthat",     "test suites"),
  c("devtools",     "package test workflow")
)
for (p in packages) {
  present <- requireNamespace(p[1], quietly = TRUE)
  version <- if (present) as.character(utils::packageVersion(p[1])) else ""
  cat(sprintf("%s\t%s\t%s\t%s\n", if (present) "ok" else "FAIL", p[1], version, p[2]))
}
EOF

  while IFS=$'\t' read -r status name version note; do
    case "${status}" in
      ok)   pass "${name}" "${version}  (${note})" ;;
      FAIL) fail "${name}" "missing -- needed for ${note}" ;;
    esac
  done <"${r_report}"
  rm -f "${r_report}"
fi

# ---------------------------------------------------------------------------
section "Python (port of the regression layer)"
# ---------------------------------------------------------------------------

if command -v python3 >/dev/null 2>&1; then
  pass "python3" "$(python3 --version 2>&1)"

  # As with R, one interpreter checks every package. `barcs` is checked last
  # and separately: it is the port itself rather than a dependency, and it
  # imports from the source tree rather than from an installed distribution.
  py_report="$(mktemp)"
  python3 - >"${py_report}" 2>/dev/null <<'EOF'
import importlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path("python").resolve()))

packages = [
    ("numpy", "arrays and linear algebra"),
    ("scipy", "distributions and the rho root-finder"),
    ("pandas", "result tables"),
    ("formulaic", "R-style formula parsing"),
    ("marimo", "the interactive notebook"),
    ("matplotlib", "notebook figures"),
    ("pytest", "the Python test suite"),
    ("barcs", "the port itself (python/barcs/)"),
]
for name, note in packages:
    try:
        module = importlib.import_module(name)
        version = getattr(module, "__version__", "")
    except Exception:
        print(f"FAIL\t{name}\t\t{note}")
    else:
        print(f"ok\t{name}\t{version}\t{note}")
EOF

  while IFS=$'\t' read -r status name version note; do
    case "${status}" in
      ok)   pass "${name}" "${version}  (${note})" ;;
      FAIL) fail "${name}" "missing -- needed for ${note}" ;;
    esac
  done <"${py_report}"
  rm -f "${py_report}"
else
  fail "python3" "not on PATH"
fi

# ---------------------------------------------------------------------------
section "Julia (CRISPulator simulations)"
# ---------------------------------------------------------------------------

if command -v julia >/dev/null 2>&1; then
  pass "julia" "$(julia --version 2>/dev/null)"
  # An installed Julia is not enough: the CRISPulator environment in `julia/`
  # has to be instantiated before `pixi run simulate` will work.
  if julia --project=julia -e 'using Crispulator' >/dev/null 2>&1; then
    pass "Crispulator.jl" "environment instantiated"
  else
    warn "Crispulator.jl" "not instantiated -- run 'pixi run julia-instantiate'"
  fi
else
  fail "julia" "not on PATH -- run 'pixi run setup-julia'"
fi

# ---------------------------------------------------------------------------
section "MAGeCK (external gene caller)"
# ---------------------------------------------------------------------------

if command -v mageck >/dev/null 2>&1; then
  mageck_version="$(mageck --version 2>&1 | tr -d '\n')"
  # Every committed MAGeCK number came from 0.5.9.5. Its gene p-value is a
  # permutation tail probability, so a different release can move results
  # without anything in this repository changing.
  if [[ "${mageck_version}" == "0.5.9.5" ]]; then
    pass "mageck" "${mageck_version}"
  else
    warn "mageck" "${mageck_version} -- committed results used 0.5.9.5"
  fi
else
  fail "mageck" "not on PATH"
fi

# ---------------------------------------------------------------------------
section "TeX (manuscript build)"
# ---------------------------------------------------------------------------

if command -v latexmk >/dev/null 2>&1 && command -v pdflatex >/dev/null 2>&1; then
  pass "latexmk" "$(pdflatex --version 2>/dev/null | head -1)"

  # cm-super is a hard requirement, not a nicety. main.tex loads microtype
  # with T1 fontenc; without scalable Type 1 Computer Modern, font expansion
  # fails and pdfTeX produces no PDF at all rather than falling back.
  if kpsewhich sfrm1000.pfb >/dev/null 2>&1; then
    pass "cm-super" "scalable Type 1 CM present"
  else
    fail "cm-super" "missing -- microtype + T1 will abort the build with no output"
  fi

  for sty in microtype.sty booktabs.sty natbib.sty amsmath.sty xcolor.sty; do
    if kpsewhich "${sty}" >/dev/null 2>&1; then
      pass "${sty%.sty}" "found"
    else
      fail "${sty%.sty}" "missing LaTeX package"
    fi
  done
else
  fail "latexmk/pdflatex" "not on PATH"
fi

# ---------------------------------------------------------------------------
section "Repository state"
# ---------------------------------------------------------------------------

# The R package under test is a submodule. An uninitialized submodule is an
# empty directory, which fails much later as a confusing R install error.
if [[ -f CB2/DESCRIPTION ]]; then
  pass "CB2 submodule" "checked out ($(git -C CB2 rev-parse --short HEAD 2>/dev/null))"
else
  fail "CB2 submodule" "empty -- run 'git submodule update --init --recursive'"
fi

# `results/` and `data/raw/` are gitignored, so a fresh clone has neither.
# That is expected, not broken: the simulations regenerate one and the
# benchmarks download the other. Reported so it is not mistaken for a fault.
if [[ -d results && -n "$(find results -mindepth 2 -type f -print -quit 2>/dev/null)" ]]; then
  pass "results/" "populated"
else
  warn "results/" "empty (gitignored) -- run 'pixi run simulate' before the FACS benchmarks"
fi

if [[ -d data/raw && -n "$(find data/raw -type f ! -name README.md -print -quit 2>/dev/null)" ]]; then
  pass "data/raw/" "populated"
else
  warn "data/raw/" "empty (gitignored) -- real-data benchmarks download their own input"
fi

# ---------------------------------------------------------------------------
printf '\n\033[1mSummary\033[0m  %d ok, %d warn, %d failed\n' \
  "${pass_count}" "${warn_count}" "${fail_count}"

if (( fail_count > 0 )); then
  printf 'Some tools are missing. Most are fixed by:  pixi run setup\n'
  exit 1
fi
printf 'Every toolchain is runnable.\n'
