#!/usr/bin/env bash
#
# Install a workspace-local TeX Live into `.tools/texlive`.
#
# This is the second dependency pixi cannot supply, for a different reason
# than Julia. conda-forge does ship `texlive-core`, but at every version it is
# only the engines: 13 MB, with an empty `texmf-dist/tex/latex`. It has no
# amsmath, no hyperref, no microtype, and its bundled `tlmgr` cannot run
# because the TeXLive Perl modules are not packaged with it. There is no
# conda-forge equivalent of Debian's texlive-latex-recommended, so a real TeX
# tree has to come from CTAN.
#
# The install starts from `scheme-infraonly` -- the installer plus tlmgr and
# nothing else -- and then adds exactly the packages main.tex loads. That
# keeps this at a few hundred MB instead of the several GB a full scheme
# costs, and it makes the dependency list auditable: every entry below is a
# `\usepackage` line in main.tex, or a documented consequence of one.
#
# Idempotent. Re-running with the tree already present only tops up packages.
#
#     pixi run setup-texlive

set -euo pipefail

ROOT="${PIXI_PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
PREFIX="${ROOT}/.tools/texlive"

# One entry per `\usepackage` in main.tex, plus three that are needed but not
# named there:
#
#   cm-super   main.tex loads microtype under T1 fontenc. Without scalable
#              Type 1 Computer Modern, microtype's font expansion fails and
#              pdfTeX exits having written no PDF at all -- it does not fall
#              back to bitmap fonts. This is the single most confusing way
#              this build can break, so it is installed explicitly.
#   ec         The T1-encoded Computer Modern fonts themselves.
#   latexmk    The build driver `pixi run manuscript` calls.
#
# tlmgr resolves each package's own dependencies, so hyperref pulls its dozen
# support packages without them being listed.
TEXLIVE_PACKAGES=(
  scheme-basic   # latex, plain, the base font set
  latexmk
  amsmath amsfonts   # amsmath, amssymb, amsthm
  booktabs
  graphics           # graphicx
  xcolor
  geometry
  hyperref
  microtype
  natbib
  cm-super
  ec
)

# `--texdir` below fixes the tree at exactly ${PREFIX}. Left to its own
# devices install-tl would append the release year (.tools/texlive/2026/),
# which would move every path in this script and in pixi.toml's PATH each
# January.
export TEXLIVE_INSTALL_TEXMFLOCAL="${PREFIX}/texmf-local"
export TEXLIVE_INSTALL_TEXMFSYSCONFIG="${PREFIX}/texmf-config"
export TEXLIVE_INSTALL_TEXMFSYSVAR="${PREFIX}/texmf-var"
# Keep per-user TeX state inside the workspace too, so a build here cannot
# pick up or write to a TeX tree elsewhere on the machine.
export TEXLIVE_INSTALL_TEXMFHOME="${PREFIX}/texmf-home"
export TEXLIVE_INSTALL_TEXMFCONFIG="${PREFIX}/texmf-config-user"
export TEXLIVE_INSTALL_TEXMFVAR="${PREFIX}/texmf-var-user"
export TEXLIVE_INSTALL_NO_WELCOME=1

find_tlmgr() {
  # TeX Live installs binaries under a platform-named directory, for example
  # bin/aarch64-linux. Rather than hardcode the mapping from `uname -m`, just
  # look for whatever the installer created.
  #
  # The `|| true` matters: on a first run `${PREFIX}` does not exist, find
  # exits nonzero, and under `set -e` that would abort the script here --
  # before the install that creates the directory has had a chance to run.
  #
  # No `-type f`: TeX Live installs bin/<platform>/tlmgr as a symlink into
  # texmf-dist/scripts, so restricting to regular files finds nothing.
  find "${PREFIX}" -path '*/bin/*' -name tlmgr 2>/dev/null | head -1 || true
}

TLMGR="$(find_tlmgr)"

if [[ -z "${TLMGR}" ]]; then
  echo "Installing TeX Live infrastructure into ${PREFIX}"
  WORKDIR="$(mktemp -d)"
  trap 'rm -rf "${WORKDIR}"' EXIT

  # mirror.ctan.org redirects to a nearby mirror and always serves the current
  # year's tlnet tree, which is what matches the current install-tl.
  curl -fsSL --retry 3 \
    -o "${WORKDIR}/install-tl-unx.tar.gz" \
    https://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz

  mkdir -p "${WORKDIR}/installer"
  tar -xzf "${WORKDIR}/install-tl-unx.tar.gz" -C "${WORKDIR}/installer" --strip-components=1

  perl "${WORKDIR}/installer/install-tl" \
    --no-interaction \
    --scheme=infraonly \
    --texdir="${PREFIX}" \
    --no-doc-install \
    --no-src-install

  TLMGR="$(find_tlmgr)"
  if [[ -z "${TLMGR}" ]]; then
    echo "setup_texlive.sh: install-tl finished but no tlmgr was produced" >&2
    exit 1
  fi
fi

echo "Installing LaTeX packages with ${TLMGR}"
# `tlmgr install` on an already-installed package prints a notice and moves
# on, which is what makes re-running this script cheap.
"${TLMGR}" install "${TEXLIVE_PACKAGES[@]}"

BINDIR="$(dirname "${TLMGR}")"
echo
echo "TeX Live ready:"
"${BINDIR}/pdflatex" --version | head -1

# Fail here rather than at the end of a manuscript build. kpsewhich resolves
# against the tree we just populated, so a missing file now means a genuinely
# missing package, not a PATH problem.
missing=()
for f in microtype.sty booktabs.sty natbib.sty amsmath.sty xcolor.sty \
         hyperref.sty geometry.sty graphicx.sty sfrm1000.pfb; do
  PATH="${BINDIR}:${PATH}" kpsewhich "${f}" >/dev/null 2>&1 || missing+=("${f}")
done

if (( ${#missing[@]} > 0 )); then
  echo "setup_texlive.sh: still missing after install: ${missing[*]}" >&2
  exit 1
fi
echo "All LaTeX packages main.tex needs are present."
