#!/usr/bin/env bash
#
# Keep the manuscript LaTeX in step across its two homes.
#
# The text has one history and two places it is served from:
#
#   origin    github.com/jeonglab-bcm/BARCS-tex  - what `overleaf/` clones from.
#                                                  Public, so anonymous clones
#                                                  and CI need no credentials,
#                                                  and GitHub indexes the prose
#                                                  for search and diffs.
#   overleaf  latex.bioinfolder.com git bridge   - where coauthors edit.
#
# This script makes both agree. It is idempotent: run it after a fresh clone to
# add the Overleaf remote, and any time afterwards to sync.
#
#   scripts/manuscript_sync.sh
#
# It fetches Overleaf, fast-forwards `overleaf/` onto it, and pushes to both
# remotes. It refuses to guess when the two histories have genuinely diverged.

set -euo pipefail

OVERLEAF_URL="https://git@latex.bioinfolder.com/git/6a6ec212d3ea56f57964b094"
GITHUB_URL="https://github.com/jeonglab-bcm/BARCS-tex.git"

if [ ! -d overleaf/.git ] && [ ! -f overleaf/.git ]; then
  echo "error: overleaf/ is not checked out. Run from the repository root:" >&2
  echo "  git submodule update --init --recursive" >&2
  exit 1
fi

git_sub() { git -C overleaf "$@"; }

# --- remotes -----------------------------------------------------------------

if ! git_sub remote get-url origin >/dev/null 2>&1; then
  git_sub remote add origin "$GITHUB_URL"
fi
if ! git_sub remote get-url overleaf >/dev/null 2>&1; then
  echo "Adding the overleaf remote."
  git_sub remote add overleaf "$OVERLEAF_URL"
fi

# --- fetch both --------------------------------------------------------------

echo "Fetching both remotes."
git_sub fetch --quiet origin main
if ! git_sub fetch --quiet overleaf main; then
  echo "error: could not fetch the Overleaf bridge. It needs credentials;" >&2
  echo "       see DEVELOPMENT.md." >&2
  exit 1
fi

local_head=$(git_sub rev-parse HEAD)
github_head=$(git_sub rev-parse origin/main)
overleaf_head=$(git_sub rev-parse overleaf/main)

echo "  local    ${local_head:0:8}"
echo "  github   ${github_head:0:8}"
echo "  overleaf ${overleaf_head:0:8}"

if [ -n "$(git_sub status --porcelain)" ]; then
  echo "error: overleaf/ has uncommitted changes. Commit or stash them first." >&2
  exit 1
fi

# --- reconcile ---------------------------------------------------------------

# Stay on a real branch, not the detached HEAD a submodule checkout leaves.
git_sub checkout --quiet main 2>/dev/null || git_sub checkout --quiet -B main "$local_head"

# Four cases, and only the last one needs a human.
if [ "$overleaf_head" = "$local_head" ]; then
  : # Overleaf has nothing new.
elif git_sub merge-base --is-ancestor "$local_head" "$overleaf_head"; then
  echo "Fast-forwarding to the Overleaf commit."
  git_sub merge --quiet --ff-only "$overleaf_head"
elif git_sub merge-base --is-ancestor "$overleaf_head" "$local_head"; then
  echo "Local commits are ahead of Overleaf; they will be pushed."
else
  echo "error: overleaf/main and the local checkout have genuinely diverged." >&2
  echo "       Both sides have commits the other lacks. Reconcile by hand in" >&2
  echo "       overleaf/ before syncing." >&2
  exit 1
fi

head=$(git_sub rev-parse HEAD)

# --- push --------------------------------------------------------------------

pushed=0
if [ "$github_head" != "$head" ]; then
  echo "Pushing to GitHub."
  git_sub push --quiet origin "HEAD:main"
  pushed=1
fi
if [ "$overleaf_head" != "$head" ]; then
  echo "Pushing to Overleaf."
  git_sub push --quiet overleaf "HEAD:main"
  pushed=1
fi

echo
if [ "$pushed" -eq 0 ]; then
  echo "Everything already agrees at ${head:0:8}."
else
  echo "Both remotes now at ${head:0:8}."
fi

# The parent repository pins a commit, so it lags until told otherwise.
pinned=$(git rev-parse HEAD:overleaf 2>/dev/null || echo none)
if [ "$pinned" != "$head" ]; then
  echo
  echo "This repository still pins ${pinned:0:8}. To advance it:"
  echo "  git add overleaf && git commit -m 'chore: update overleaf submodule'"
fi
