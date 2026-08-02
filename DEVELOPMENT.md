# Developing the manuscript and the BARCS package

This repository has two coordinated histories:

- `jeonglab-bcm/BARCS-manuscript` (this repository) contains the manuscript,
  benchmark code, compact result summaries, and figures.
- `jeonglab-bcm/BARCS` contains the R package. It is checked out here as the
  `BARCS/` Git submodule.

A third submodule, `CB2/`, pins the older CB2 package. It is **not** a
dependency of BARCS. It is kept for two reasons only:
`examples/barcs_input_output_examples.R` runs original CB2 as a benchmark
baseline, and two scripts load the Sanson screen from `CB2/data/`.

## Clone both repositories

```sh
git clone --recurse-submodules https://github.com/jeonglab-bcm/BARCS-manuscript.git
cd BARCS-manuscript
```

For a clone created without `--recurse-submodules`:

```sh
git submodule update --init --recursive
```

## How scripts find the package

Every analysis script starts with

```r
source(file.path("R", "load_barcs.R"))
```

run from the repository root. That loader prefers an installed BARCS and
otherwise loads the `BARCS/` submodule in place with `pkgload::load_all()`, so
a fresh clone works with no install step. Installing is faster, because the
compiled kernels are then built once rather than once per session:

```sh
R CMD INSTALL BARCS
```

The loader also enforces a minimum package version, so a stale submodule fails
with a clear message rather than a missing-function error.

## Change the BARCS package

```sh
cd BARCS
git switch main

# Edit, document, and test the package.
Rscript -e 'devtools::document()'
Rscript -e 'devtools::test()'

git add <changed-package-files>
git commit -m "Describe the package change"
git push
cd ..
```

Bumping `Version:` in `BARCS/DESCRIPTION` on `main` makes the package's release
workflow tag that version and publish a source tarball. Editing DESCRIPTION
without changing the version does nothing, so the workflow is safe to re-run.

This repository still points at the previous package commit until its submodule
pointer is updated:

```sh
git add BARCS
git commit -m "chore: update BARCS submodule"
git push
```

This two-commit workflow is intentional: package code stays reviewable in its
own repository, while each analysis revision records the exact package commit
it used.

## Pull coordinated updates

```sh
git pull
git submodule update --init --recursive
```

To advance a submodule to the newest commit on its configured branch:

```sh
git submodule update --remote BARCS
git add BARCS
git commit -m "chore: update BARCS submodule"
```

## Verify before publishing

```sh
Rscript -e 'devtools::test("BARCS")'
R CMD build BARCS
R CMD check --no-manual BARCS_*.tar.gz
Rscript examples/barcs_quickstart.R
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The package's own test suite lives in `BARCS/tests/testthat/`. It replaced the
former `tests/run_tests.R` in this repository when BARCS was split out, so
there is one copy of both the implementation and its tests.
