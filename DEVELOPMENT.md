# Developing the parent project and CB2

This repository has two coordinated histories:

- `jeonglab-bcm/BARCS` contains the manuscript, benchmark code, compact
  result summaries, and figures.
- `jeonglab-bcm/CB2` contains the R package. It is checked out here as the
  `CB2/` Git submodule.

## Clone both repositories

```sh
git clone --recurse-submodules https://github.com/jeonglab-bcm/BARCS.git
cd BARCS
```

For a clone created without `--recurse-submodules`:

```sh
git submodule update --init --recursive
```

## Change the CB2 package

The parent project currently follows the `feature/barcs` package branch.

```sh
cd CB2
git switch feature/barcs

# Edit and test the package.
Rscript -e 'devtools::test()'

git add <changed-package-files>
git commit -m "Describe the package change"
git push
cd ..
```

The parent repository still points to the previous CB2 commit until its
submodule pointer is updated:

```sh
git add CB2
git commit -m "chore: update CB2 submodule"
git push
```

This two-commit workflow is intentional: package code remains reviewable in
CB2, while each analysis revision records the exact package commit it used.

## Pull coordinated updates

```sh
git pull
git submodule update --init --recursive
```

To advance the submodule to the newest commit on its configured branch:

```sh
git submodule update --remote CB2
git add CB2
git commit -m "chore: update CB2 submodule"
```

## Verify before publishing

```sh
Rscript tests/run_tests.R
Rscript -e 'devtools::test("CB2")'
R CMD build CB2
R CMD check --no-manual CB2_*.tar.gz
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```
