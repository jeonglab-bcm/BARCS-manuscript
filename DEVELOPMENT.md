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

## Set up the toolchain

The analysis spans four toolchains, and reproducing a benchmark means having
all four at the right versions. They are declared in
[`pixi.toml`](pixi.toml) so that is one command rather than an afternoon:

```sh
curl -fsSL https://pixi.sh/install.sh | bash   # once, if pixi is missing
pixi run setup
pixi run doctor
```

`pixi run doctor` prints one line per capability — R and each package, Julia
and the CRISPulator environment, MAGeCK and its version, every LaTeX package
`main.tex` loads, and whether the submodule is checked out. Run it first when
something fails for a reason that looks environmental.

Work either through tasks (`pixi run test`, `pixi task list` for the rest) or
inside `pixi shell`, where the whole toolchain is on `PATH` and plain
`Rscript`, `julia`, `mageck`, and `latexmk` all work.

### What pixi installs, and what it cannot

Most of it is conda packages pinned by `pixi.lock`: R, the CRAN and
Bioconductor packages, the compilers, and MAGeCK 0.5.9.5 from bioconda. Three
things need separate handling, each for a stated reason:

| Component | Why it is not a conda dependency | Installed by |
|---|---|---|
| Julia 1.12.6 | conda-forge builds Julia for linux-64 only; this project also runs on aarch64. Fetched from the official binaries against a pinned checksum. | `scripts/setup_julia.sh` |
| TeX Live | conda-forge's `texlive-core` is 13 MB of engines with an empty `texmf-dist/tex/latex` — no amsmath, no hyperref, no microtype — and its bundled `tlmgr` cannot run. Installed from CTAN, starting at `scheme-infraonly` and adding exactly what `main.tex` loads. | `scripts/setup_texlive.sh` |
| CRISPhieRmix, simCRISPR | Neither is on CRAN, Bioconductor, or any conda channel. Pinned by commit instead. | `scripts/setup_r_github.R` |

Both external toolchains land in the gitignored `.tools/`, and every R package
goes into the pixi environment's own library, so nothing here touches a
user-level R library or `~/.julia`.

> **The manuscript build needs `cm-super`.** `main.tex` loads `microtype`
> under `T1` `fontenc`. Without scalable Type 1 Computer Modern, microtype's
> font expansion fails and pdfTeX exits having written **no PDF at all** — it
> does not fall back to bitmap fonts, and the failure does not name the
> missing font. `setup_texlive.sh` installs it and `pixi run doctor` checks
> for it. On a system TeX Live, `texlive-latex-recommended`, `-extra`,
> `-fonts-recommended`, and `-science` are not sufficient on their own.

### Version pins that are not cosmetic

Three pins exist because moving them changes committed numbers, not just
build metadata:

- **MAGeCK 0.5.9.5.** Its gene p-value is a permutation tail probability, so a
  different release moves results. `R/mageck.R` warns when it sees another
  version.
- **Julia 1.12.6.** CRISPulator draws screens from Julia's global RNG, whose
  stream is only guaranteed within a release series. A different Julia gives
  different screens.
- **CRISPhieRmix and simCRISPR commits.** Both feed gene calls or simulated
  data into committed metrics.

## Change the CB2 package

The submodule tracks the `hj/manuscript-audit` branch, as recorded in
`.gitmodules`.

```sh
cd CB2
git switch hj/manuscript-audit

# Edit and test the package.
pixi run --manifest-path .. test-package    # or: Rscript -e 'devtools::test()'

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

After changing the package, reinstall it into the pixi environment so the
benchmarks pick the change up:

```sh
pixi run setup-cb2
```

Note that `R/bbreg.R` in this repository is a standalone copy of the
regression layer, used so a benchmark script can be read and run without the
package installed. A change to the package's regression code has to be
mirrored there, and `pixi run test` covers that copy while
`pixi run test-package` covers the package.

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
pixi run doctor          # every toolchain present
pixi run test            # analysis-side regression tests
pixi run test-package    # the package's own suite
pixi run check-package   # R CMD build + R CMD check
pixi run manuscript      # latexmk; must finish with no undefined references
```

Long-form equivalents, if you are not using pixi:

```sh
Rscript tests/run_tests.R
Rscript -e 'devtools::test("CB2")'
R CMD build CB2
R CMD check --no-manual CB2_*.tar.gz
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```
