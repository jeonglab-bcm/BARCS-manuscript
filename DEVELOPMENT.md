# Developing the manuscript and the BARCS package

This repository coordinates three histories of its own plus one legacy pin:

- `jeonglab-bcm/BARCS-manuscript` (this repository) holds the benchmark code,
  compact result summaries, and generated figures.
- `jeonglab-bcm/BARCS` holds the R package, checked out as the `BARCS/`
  submodule.
- `jeonglab-bcm/BARCS-tex` holds the manuscript LaTeX, checked out as the
  `overleaf/` submodule. It is a public mirror of the Overleaf project where
  the writing actually happens. See "Change the manuscript text" below.

`CB2/` pins the older CB2 package. It is **not** a dependency of BARCS. It is
kept for two reasons only: `examples/barcs_input_output_examples.R` runs
original CB2 as a benchmark baseline, and two scripts load the Sanson screen
from `CB2/data/`.

Nothing is duplicated between them: each artifact has exactly one home, and
this repository pins the commit of each that produced a given revision.

## Clone the repository and its submodules

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

## Change the manuscript text

The manuscript has one history served from two places:

| Remote | Where | Why |
|---|---|---|
| `origin` | `github.com/jeonglab-bcm/BARCS-tex` | What `overleaf/` clones from. Public, so anonymous clones and CI need no credentials, and **GitHub indexes the prose for search and diffs**. |
| `overleaf` | `latex.bioinfolder.com` Git bridge | Where coauthors edit. Authenticated. |

`scripts/manuscript_sync.sh` keeps them in step. It is idempotent: run it after
a fresh clone to add the Overleaf remote, and any time afterwards to sync.

```sh
scripts/manuscript_sync.sh          # fetch Overleaf, fast-forward, push to both
$EDITOR overleaf/sections/3_barcs-results.tex
git -C overleaf commit -am "Revise the results section"
scripts/manuscript_sync.sh          # push the edit to both remotes

git add overleaf
git commit -m "chore: update overleaf submodule"
```

The script refuses to guess: it stops if `overleaf/` is dirty, or if the two
remotes have genuinely diverged with commits on both sides.

Note the submodule URL changed when the mirror was introduced. An existing
clone needs one command to notice:

```sh
git submodule sync -- overleaf
```

Build the PDF from the submodule; `-cd` makes latexmk run in `overleaf/` so the
relative `\input` and `figures/` paths resolve:

```sh
latexmk -pdf -cd overleaf/main.tex
```

### Publishing a regenerated figure

Analysis scripts write every figure to `figures/` in this repository, but only
the three that `\includegraphics` names belong in the Overleaf project. After
regenerating one:

```sh
Rscript examples/manuscript_liang_figure.R
scripts/publish_figures.sh
git -C overleaf add figures && git -C overleaf commit -m "Update figures"
git -C overleaf push
git add overleaf && git commit -m "chore: update overleaf submodule"
```

`scripts/publish_figures.sh` copies only those three, reports what changed, and
fails if a cited figure has not been built. Everything else in `figures/` is a
diagnostic plot and is gitignored. Follow it with
`scripts/manuscript_sync.sh` to push the new figure to both remotes.

### Continuous PDF builds

`.github/workflows/manuscript.yaml` compiles the manuscript and fails if the
LaTeX does not build, attaching the PDF as an artifact and the LaTeX log when a
build fails. It runs when the `overleaf/` submodule pointer moves on `main` or
in a pull request, and weekly on a schedule.

The schedule matters because Overleaf edits do not touch this repository at
all: until someone bumps the submodule pointer, this repository still pins the
old commit. The weekly run builds the newest commit on the Overleaf side
instead of the pinned one, so a coauthor's broken LaTeX surfaces without
waiting for a pointer bump. `workflow_dispatch` can do the same on demand with
the `use_latest` input.

No secret is required, because the submodule clones from the public mirror.

The `OVERLEAF_GIT_TOKEN` secret is optional and only affects the scheduled run:
with it, that run builds straight from the Overleaf bridge and so catches a
broken edit before anyone has mirrored it. Without it, the scheduled run builds
the mirror's newest commit instead, which is still ahead of the pinned one.

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
latexmk -pdf -cd -interaction=nonstopmode -halt-on-error overleaf/main.tex
```

The package's own test suite lives in `BARCS/tests/testthat/`. It replaced the
former `tests/run_tests.R` in this repository when BARCS was split out, so
there is one copy of both the implementation and its tests.
