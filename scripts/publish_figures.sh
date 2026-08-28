#!/usr/bin/env bash
#
# Copy the figures the manuscript cites into the Overleaf submodule.
#
# The manuscript lives in `overleaf/`, which is the Overleaf project's Git
# bridge. Analysis scripts write every figure into `figures/`, but only the
# three below are included by \includegraphics, so only those belong in the
# Overleaf project.
#
# Run from the repository root:
#
#   scripts/publish_figures.sh
#
# Then commit and push inside `overleaf/` to make the change appear in
# Overleaf, and update this repository's submodule pointer.

set -euo pipefail

CITED=(
  crispulator_facs_moi_10k_f1_by_fdr.pdf
  liang_longitudinal_volcano_trajectories.pdf
  simcrispr_interaction.pdf
  liu_tcell_publication_comparison.pdf
)

if [ ! -d overleaf ]; then
  echo "error: overleaf/ is missing. Run from the repository root, after:" >&2
  echo "  git submodule update --init --recursive" >&2
  exit 1
fi

mkdir -p overleaf/figures

missing=0
copied=0
for figure in "${CITED[@]}"; do
  source_path="figures/${figure}"
  target_path="overleaf/figures/${figure}"

  if [ ! -f "$source_path" ]; then
    echo "missing   ${source_path} (regenerate it; see figures/README.md)" >&2
    missing=$((missing + 1))
    continue
  fi

  if [ -f "$target_path" ] && cmp -s "$source_path" "$target_path"; then
    echo "unchanged ${figure}"
    continue
  fi

  cp "$source_path" "$target_path"
  echo "copied    ${figure}"
  copied=$((copied + 1))
done

if [ "$missing" -gt 0 ]; then
  echo >&2
  echo "error: ${missing} cited figure(s) are not built; nothing was published for them." >&2
  exit 1
fi

echo
if [ "$copied" -eq 0 ]; then
  echo "Overleaf already has the current figures."
else
  echo "Copied ${copied} figure(s). To publish them:"
  echo "  git -C overleaf add figures && git -C overleaf commit -m 'Update figures' && git -C overleaf push"
  echo "  git add overleaf && git commit -m 'chore: update overleaf submodule'"
fi
